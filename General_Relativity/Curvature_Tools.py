"""Curvature_Tools.py - the tensor calculus shared by every metric in this project.

Given a metric g_{mu nu} and a list of coordinates, this module works out

    * the Christoffel symbols   Gamma^l_{ij}
    * the Ricci tensor          R_{ij}
    * the Ricci scalar          R
    * the Riemann tensor and the Kretschmann scalar  K = R_{abcd} R^{abcd}

Conventions are the usual ones (Misner-Thorne-Wheeler): signature (-, +, +, +),
c = G = 1, and

    Gamma^l_{ij}   = 1/2 g^{lk} (d_i g_{kj} + d_j g_{ki} - d_k g_{ij})
    R^r_{s m n}    = d_m Gamma^r_{ns} - d_n Gamma^r_{ms}
                     + Gamma^r_{ml} Gamma^l_{ns} - Gamma^r_{nl} Gamma^l_{ms}
    R_{s n}        = R^r_{s r n}

Nothing here depends on four dimensions, so the same functions work for any
metric you care to type in.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Optional

import sympy as sp


@dataclass
class Spacetime:
    """A metric together with the coordinates it is written in."""

    name: str
    metric: sp.Matrix
    coords: list
    description: str = ""
    # If you know the inverse metric in closed form, pass it here. It saves a slow
    # symbolic matrix inversion (this matters for Kerr).
    inverse: Optional[sp.Matrix] = None
    # Function used to tidy up every result. sp.simplify is thorough but slow.
    simplify: Callable = sp.simplify

    @property
    def dim(self) -> int:
        return len(self.coords)


@dataclass
class Curvature:
    """Everything computed for one spacetime."""

    spacetime: Spacetime
    christoffel: list                      # christoffel[l][i][j] = Gamma^l_{ij}
    ricci: sp.Matrix
    ricci_scalar: sp.Expr
    kretschmann: Optional[sp.Expr] = None
    timings: dict = field(default_factory=dict)

    @property
    def is_ricci_flat(self) -> bool:
        return all(entry == 0 for entry in self.ricci)


# ---------------------------------------------------------------------------
# The calculations
# ---------------------------------------------------------------------------

def inverse_metric(spacetime: Spacetime) -> sp.Matrix:
    """The inverse metric g^{mu nu}, using the supplied closed form if there is one."""
    if spacetime.inverse is not None:
        return spacetime.inverse
    return spacetime.metric.inv().applyfunc(spacetime.simplify)


def christoffel_symbols(spacetime: Spacetime, g_inv: Optional[sp.Matrix] = None) -> list:
    """Christoffel symbols of the second kind, as christoffel[l][i][j] = Gamma^l_{ij}."""
    g, x, n = spacetime.metric, spacetime.coords, spacetime.dim
    g_inv = inverse_metric(spacetime) if g_inv is None else g_inv
    simplify = spacetime.simplify

    # dg[k][i][j] = d g_{ij} / d x^k, worked out once instead of n^3 times
    dg = [[[sp.diff(g[i, j], x[k]) for j in range(n)] for i in range(n)] for k in range(n)]

    gamma = [[[sp.S.Zero] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for i in range(n):
            for j in range(i, n):                     # Gamma^l_{ij} = Gamma^l_{ji}
                total = sum(
                    g_inv[l, k] * (dg[i][k][j] + dg[j][k][i] - dg[k][i][j])
                    for k in range(n) if g_inv[l, k] != 0
                )
                value = simplify(sp.Rational(1, 2) * total)
                gamma[l][i][j] = gamma[l][j][i] = value
    return gamma


def ricci_tensor(spacetime: Spacetime, gamma: list) -> sp.Matrix:
    """R_{sn} = d_r Gamma^r_{ns} - d_n Gamma^r_{rs}
                + Gamma^r_{rl} Gamma^l_{ns} - Gamma^r_{nl} Gamma^l_{rs}

    The last two products are different sums (they have different index
    structure); the Ricci tensor needs both.
    """
    x, n, simplify = spacetime.coords, spacetime.dim, spacetime.simplify

    # contracted symbols Gamma^r_{rl}: handy and used a lot
    trace = [sum(gamma[r][r][l] for r in range(n)) for l in range(n)]

    ricci = sp.zeros(n, n)
    for s in range(n):
        for nu in range(s, n):                        # the Ricci tensor is symmetric
            value = (
                sum(sp.diff(gamma[r][nu][s], x[r]) for r in range(n))
                - sp.diff(trace[s], x[nu])
                + sum(trace[l] * gamma[l][nu][s] for l in range(n))
                - sum(gamma[r][nu][l] * gamma[l][r][s] for r in range(n) for l in range(n))
            )
            ricci[s, nu] = ricci[nu, s] = simplify(value)
    return ricci


def ricci_scalar(spacetime: Spacetime, ricci: sp.Matrix, g_inv: Optional[sp.Matrix] = None) -> sp.Expr:
    """R = g^{mu nu} R_{mu nu}"""
    g_inv = inverse_metric(spacetime) if g_inv is None else g_inv
    n = spacetime.dim
    return spacetime.simplify(sum(g_inv[i, j] * ricci[i, j] for i in range(n) for j in range(n)))


def riemann_lowered(spacetime: Spacetime, gamma: list) -> dict:
    """Every non-zero component R_{abcd} of the Riemann tensor, as {(a, b, c, d): expression}.

    Components related by the antisymmetry in the last two indices are both
    stored, so a contraction can just loop over the dictionary.
    """
    g, x, n, simplify = spacetime.metric, spacetime.coords, spacetime.dim, spacetime.simplify

    def mixed(r, s, m, nu):                           # R^r_{s m nu}
        return (sp.diff(gamma[r][nu][s], x[m]) - sp.diff(gamma[r][m][s], x[nu])
                + sum(gamma[r][m][l] * gamma[l][nu][s] - gamma[r][nu][l] * gamma[l][m][s]
                      for l in range(n)))

    components = {}
    for s in range(n):
        for m in range(n):
            for nu in range(m + 1, n):
                column = [mixed(r, s, m, nu) for r in range(n)]
                for a in range(n):                    # lower the first index with g_{a r}
                    value = simplify(sum(g[a, r] * column[r] for r in range(n) if g[a, r] != 0))
                    if value != 0:
                        components[(a, s, m, nu)] = value
                        components[(a, s, nu, m)] = -value
    return components


def kretschmann_scalar(spacetime: Spacetime, gamma: list, g_inv: Optional[sp.Matrix] = None) -> sp.Expr:
    """K = R_{abcd} R^{abcd}, the standard coordinate-independent curvature invariant."""
    g_inv = inverse_metric(spacetime) if g_inv is None else g_inv
    n = spacetime.dim
    low = riemann_lowered(spacetime, gamma)

    total = sp.S.Zero
    for (a, b, c, d), left in low.items():
        for (p, q, r, s), right in low.items():
            weight = g_inv[a, p] * g_inv[b, q] * g_inv[c, r] * g_inv[d, s]
            if weight != 0:
                total += weight * left * right
    return spacetime.simplify(total)


def compute_curvature(spacetime: Spacetime, with_kretschmann: bool = False,
                      progress: Optional[Callable[[str], None]] = None) -> Curvature:
    """Run the whole pipeline. `progress` gets a short message before each stage."""
    say = progress or (lambda message: None)
    timings = {}

    def timed(label, func):
        say(label)
        start = time.perf_counter()
        result = func()
        timings[label] = time.perf_counter() - start
        return result

    g_inv = timed("inverse metric", lambda: inverse_metric(spacetime))
    gamma = timed("Christoffel symbols", lambda: christoffel_symbols(spacetime, g_inv))
    ricci = timed("Ricci tensor", lambda: ricci_tensor(spacetime, gamma))
    scalar = timed("Ricci scalar", lambda: ricci_scalar(spacetime, ricci, g_inv))
    kretschmann = None
    if with_kretschmann:
        kretschmann = timed("Kretschmann scalar", lambda: kretschmann_scalar(spacetime, gamma, g_inv))
    return Curvature(spacetime, gamma, ricci, scalar, kretschmann, timings)


# ---------------------------------------------------------------------------
# Presentation
# ---------------------------------------------------------------------------

def nonzero_christoffel(curvature: Curvature) -> list:
    """[(l, i, j, expression), ...] for i <= j, skipping zeros and the mirror images."""
    n = curvature.spacetime.dim
    return [(l, i, j, curvature.christoffel[l][i][j])
            for l in range(n) for i in range(n) for j in range(i, n)
            if curvature.christoffel[l][i][j] != 0]


def format_matrix(matrix: sp.Matrix, latex: bool = False) -> str:
    """One row per line, so wide expressions stay readable."""
    if latex:
        return sp.latex(matrix)
    lines = []
    for i in range(matrix.rows):
        lines.append("  row %d: %s" % (i, ",   ".join(str(matrix[i, j]) for j in range(matrix.cols))))
    return "\n".join(lines)


def print_report(curvature: Curvature, latex: bool = False, show_christoffel: bool = True) -> None:
    """Print the metric, Christoffel symbols, Ricci tensor and scalar."""
    st = curvature.spacetime
    names = [str(c) for c in st.coords]
    show = (lambda e: sp.latex(e)) if latex else str

    print(f"\n{'=' * 70}\n{st.name}\n{'=' * 70}")
    if st.description:
        print(st.description)
    print(f"Coordinates: ({', '.join(names)})")

    print("\nMetric g_{mu nu}:")
    for i in range(st.dim):
        for j in range(i, st.dim):
            if st.metric[i, j] != 0:
                print(f"  g[{names[i]},{names[j]}] = {show(st.metric[i, j])}")

    if show_christoffel:
        print("\nChristoffel symbols (non-zero, Gamma^l_{ij} = Gamma^l_{ji} so only i <= j is listed):")
        for l, i, j, expr in nonzero_christoffel(curvature):
            print(f"  Gamma^{names[l]}_({names[i]},{names[j]}) = {show(expr)}")

    print("\nRicci tensor R_{mu nu} (non-zero components):")
    any_ricci = False
    for i in range(st.dim):
        for j in range(i, st.dim):
            if curvature.ricci[i, j] != 0:
                any_ricci = True
                print(f"  R[{names[i]},{names[j]}] = {show(curvature.ricci[i, j])}")
    if not any_ricci:
        print("  all components are zero (Ricci-flat: a vacuum solution)")

    print(f"\nRicci scalar R = {show(curvature.ricci_scalar)}")
    if curvature.kretschmann is not None:
        print(f"Kretschmann scalar K = {show(curvature.kretschmann)}")

    total = sum(curvature.timings.values())
    print(f"\n(computed in {total:.1f} s)")
