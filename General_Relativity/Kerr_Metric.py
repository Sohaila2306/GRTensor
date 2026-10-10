"""Kerr_Metric.py - the spacetime around a spinning black hole.

Roy Kerr found this exact solution in 1963. It is written here in
Boyer-Lindquist coordinates (1967), with G = c = 1:

    rho^2 = r^2 + a^2 cos^2(theta)          Delta = r^2 - 2 M r + a^2

    ds^2 = -(1 - 2 M r / rho^2) dt^2 - (4 M a r sin^2(theta) / rho^2) dt dphi
           + (rho^2 / Delta) dr^2 + rho^2 dtheta^2
           + (r^2 + a^2 + 2 M r a^2 sin^2(theta) / rho^2) sin^2(theta) dphi^2

M is the mass and a = J/M is the spin per unit mass. Setting a = 0 gives Schwarzschild.

Kerr is much heavier than the other two metrics because it has an off-diagonal
term and everything depends on both r and theta. SymPy's general-purpose
simplify() takes many minutes on it, so this file uses two shortcuts:

  * the inverse metric is written down in closed form instead of being computed, and
  * sin(theta) and cos(theta) are treated as ordinary algebraic variables S and C
    that obey S^2 + C^2 = 1, which lets SymPy use fast polynomial arithmetic.

With those, the Christoffel symbols take about a second and the Ricci tensor
(which comes out as exactly zero) about ten.

Run this file directly to see the curvature:    python Kerr_Metric.py
"""

import sympy as sp

from Curvature_Tools import Spacetime, compute_curvature, print_report

t, r, theta, phi = sp.symbols("t r theta phi", real=True)
M = sp.symbols("M", positive=True)          # mass
a = sp.symbols("a", real=True)              # spin per unit mass, J / M

coords = [t, r, theta, phi]

sin_t, cos_t = sp.sin(theta), sp.cos(theta)
rho2 = r**2 + a**2 * cos_t**2
delta = r**2 - 2 * M * r + a**2

g = sp.zeros(4, 4)
g[0, 0] = -(1 - 2 * M * r / rho2)
g[0, 3] = g[3, 0] = -2 * M * a * r * sin_t**2 / rho2
g[1, 1] = rho2 / delta
g[2, 2] = rho2
g[3, 3] = (r**2 + a**2 + 2 * M * r * a**2 * sin_t**2 / rho2) * sin_t**2

# Inverse metric in closed form (checked against g in the tests).
g_inv = sp.zeros(4, 4)
g_inv[0, 0] = -((r**2 + a**2) ** 2 - a**2 * delta * sin_t**2) / (rho2 * delta)
g_inv[0, 3] = g_inv[3, 0] = -2 * M * a * r / (rho2 * delta)
g_inv[1, 1] = delta / rho2
g_inv[2, 2] = 1 / rho2
g_inv[3, 3] = (delta - a**2 * sin_t**2) / (rho2 * delta * sin_t**2)


# --- fast simplification -----------------------------------------------------

_S, _C = sp.symbols("S C", real=True)       # stand-ins for sin(theta) and cos(theta)


def _reduce_trig(poly):
    """Rewrite a polynomial in S and C using S^2 = 1 - C^2, so S never appears squared."""
    if poly.is_number:
        return poly
    return sp.reduced(sp.expand(poly), [_S**2 + _C**2 - 1], _S, _C, order="lex")[1]


def fast_simplify(expr):
    """Much quicker than sp.simplify for expressions rational in r, sin(theta), cos(theta)."""
    e = sp.sympify(expr).subs({sin_t: _S, cos_t: _C})
    e = sp.cancel(sp.together(e))
    numerator, denominator = sp.fraction(e)
    if numerator == 0:
        return sp.S.Zero
    e = sp.cancel(_reduce_trig(numerator) / _reduce_trig(denominator))
    return e.subs({_S: sin_t, _C: cos_t})


def spacetime() -> Spacetime:
    return Spacetime(
        name="Kerr metric (Boyer-Lindquist coordinates)",
        metric=g,
        coords=coords,
        description="Vacuum around a rotating black hole with mass M and spin a = J/M (G = c = 1).",
        inverse=g_inv,
        simplify=fast_simplify,
    )


if __name__ == "__main__":
    # No Kretschmann scalar here: symbolically it takes too long. See the README for its closed form.
    print_report(compute_curvature(spacetime()))
