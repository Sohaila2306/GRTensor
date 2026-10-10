"""Tests for the curvature code. Run with:  python -m pytest   (or just: python test_general_relativity.py)

Every expected value below is a standard textbook result, so these tests check
the physics and not just that the code runs. The Kerr test takes about ten seconds.
"""

import io
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

import sympy as sp

import FLRW_Metric
import GeneralRelativity_Project as project
import Kerr_Metric
import Schwarzschild_Metric
from Curvature_Tools import Spacetime, compute_curvature


def is_zero(expr):
    return sp.simplify(expr) == 0


# --- the bug in the original code ---------------------------------------------

def test_sphere_ricci_tensor_needs_the_quadratic_christoffel_terms():
    """A 2-sphere of radius a has R_ij = g_ij / a^2 and R = 2/a^2.

    The original Ricci code failed this: its two Gamma*Gamma terms were the same sum
    with the indices renamed, so they cancelled and half of the formula was missing.
    (Its scalar happened to come out right for the sphere because the errors cancel in
    the trace, so the test has to look at the tensor components.)
    """
    th, ph = sp.symbols("theta phi", real=True)
    a = sp.symbols("a", positive=True)
    g = sp.diag(a**2, a**2 * sp.sin(th) ** 2)
    result = compute_curvature(Spacetime("sphere", g, [th, ph]), with_kretschmann=True)
    assert (result.ricci - g / a**2).applyfunc(sp.simplify) == sp.zeros(2, 2)
    assert is_zero(result.ricci_scalar - 2 / a**2)
    assert is_zero(result.kretschmann - 4 / a**4)


def test_flat_space_in_spherical_coordinates_has_no_curvature():
    t, r, th, ph = sp.symbols("t r theta phi", real=True)
    flat = Spacetime("flat", sp.diag(-1, 1, r**2, r**2 * sp.sin(th) ** 2), [t, r, th, ph])
    result = compute_curvature(flat, with_kretschmann=True)
    assert result.is_ricci_flat and result.ricci_scalar == 0 and result.kretschmann == 0


def test_no_floating_point_numbers_leak_into_exact_results():
    result = compute_curvature(Schwarzschild_Metric.spacetime())
    for l in range(4):
        for i in range(4):
            for j in range(4):
                assert not result.christoffel[l][i][j].has(sp.Float)


# --- Schwarzschild ---------------------------------------------------------------

def test_schwarzschild_christoffel_symbols():
    S = Schwarzschild_Metric
    result = compute_curvature(S.spacetime())
    G = result.christoffel
    M, r, th = S.M, S.r, S.theta
    assert is_zero(G[0][0][1] - M / (r * (r - 2 * M)))                # Gamma^t_tr
    assert is_zero(G[1][0][0] - M * (r - 2 * M) / r**3)               # Gamma^r_tt
    assert is_zero(G[1][1][1] + M / (r * (r - 2 * M)))                # Gamma^r_rr
    assert is_zero(G[1][2][2] + (r - 2 * M))                          # Gamma^r_thth
    assert is_zero(G[1][3][3] + (r - 2 * M) * sp.sin(th) ** 2)        # Gamma^r_phph
    assert is_zero(G[2][1][2] - 1 / r) and is_zero(G[3][1][3] - 1 / r)
    assert is_zero(G[2][3][3] + sp.sin(th) * sp.cos(th))
    assert is_zero(G[3][2][3] - sp.cos(th) / sp.sin(th))


def test_schwarzschild_is_a_vacuum_solution():
    result = compute_curvature(Schwarzschild_Metric.spacetime(), with_kretschmann=True)
    assert result.is_ricci_flat
    assert result.ricci_scalar == 0


def test_schwarzschild_kretschmann_scalar():
    S = Schwarzschild_Metric
    result = compute_curvature(S.spacetime(), with_kretschmann=True)
    assert is_zero(result.kretschmann - 48 * S.M**2 / S.r**6)
    # finite at the horizon r = 2M (a coordinate singularity), infinite at r = 0 (a real one)
    assert sp.simplify(result.kretschmann.subs(S.r, 2 * S.M)) == sp.Rational(3, 4) / S.M**4


def test_christoffel_symbols_are_symmetric_in_their_lower_indices():
    result = compute_curvature(FLRW_Metric.spacetime())
    for l in range(4):
        for i in range(4):
            for j in range(4):
                assert result.christoffel[l][i][j] == result.christoffel[l][j][i]
    assert result.ricci == result.ricci.T


# --- writing your own metric (the README example) ----------------------------------

def test_reissner_nordstrom_charged_black_hole():
    t, r, theta, phi = sp.symbols("t r theta phi", real=True)
    M, Q = sp.symbols("M Q", positive=True)
    f = 1 - 2 * M / r + Q**2 / r**2
    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta) ** 2)
    result = compute_curvature(Spacetime("Reissner-Nordstrom", g, [t, r, theta, phi]), with_kretschmann=True)
    assert result.ricci_scalar == 0                      # electromagnetic fields are traceless
    assert not result.is_ricci_flat                      # but the field carries energy, so R_mu_nu != 0
    assert is_zero(result.kretschmann - 8 * (6 * M**2 * r**2 - 12 * M * Q**2 * r + 7 * Q**4) / r**8)


# --- FLRW ------------------------------------------------------------------------

def test_flrw_matches_the_friedmann_equations_ingredients():
    F = FLRW_Metric
    result = compute_curvature(F.spacetime(), with_kretschmann=True)
    a, k, t = F.a, F.k, F.t
    ad, add = sp.diff(a, t), sp.diff(a, t, 2)
    assert is_zero(result.ricci[0, 0] + 3 * add / a)                                   # R_tt
    assert is_zero(result.ricci_scalar - 6 * (add / a + ad**2 / a**2 + k / a**2))      # R
    assert is_zero(result.kretschmann - 12 * ((add / a) ** 2 + ((ad**2 + k) / a**2) ** 2))


def test_de_sitter_space():
    F = FLRW_Metric
    H = sp.symbols("H", positive=True)
    result = compute_curvature(F.spacetime(0, sp.exp(H * F.t)), with_kretschmann=True)
    assert is_zero(result.ricci_scalar - 12 * H**2)
    assert is_zero(result.kretschmann - 24 * H**4)
    # de Sitter is an "Einstein space": R_mu_nu = 3 H^2 g_mu_nu
    difference = (result.ricci - 3 * H**2 * F.metric(0, sp.exp(H * F.t))).applyfunc(sp.simplify)
    assert difference == sp.zeros(4, 4)


def test_flrw_with_fixed_curvature_matches_the_general_result():
    F = FLRW_Metric
    general = compute_curvature(F.spacetime())
    closed = compute_curvature(F.spacetime(k_value=1))
    assert is_zero(closed.ricci_scalar - general.ricci_scalar.subs(F.k, 1))


def test_matter_dominated_universe_has_R_equal_to_8_pi_rho():
    """For dust a(t) ~ t^(2/3) and rho = 1/(6 pi t^2), so R = 8 pi rho = 4 / (3 t^2) (G = c = 1)."""
    F = FLRW_Metric
    result = compute_curvature(F.spacetime(0, F.t ** sp.Rational(2, 3)))
    assert is_zero(result.ricci_scalar - sp.Rational(4, 3) / F.t**2)


# --- Kerr ------------------------------------------------------------------------

def test_kerr_inverse_metric_really_is_the_inverse():
    product = (Kerr_Metric.g * Kerr_Metric.g_inv - sp.eye(4)).applyfunc(Kerr_Metric.fast_simplify)
    assert product == sp.zeros(4, 4)


def test_kerr_reduces_to_schwarzschild_when_a_is_zero():
    kerr = Kerr_Metric.g.subs(Kerr_Metric.a, 0).applyfunc(sp.simplify)
    expected = sp.diag(*[Schwarzschild_Metric.g[i, i] for i in range(4)]).subs(
        {Schwarzschild_Metric.M: Kerr_Metric.M, Schwarzschild_Metric.r: Kerr_Metric.r,
         Schwarzschild_Metric.theta: Kerr_Metric.theta})
    assert (kerr - expected).applyfunc(sp.simplify) == sp.zeros(4, 4)


def test_kerr_is_a_vacuum_solution_and_has_the_known_gamma_r_tt():
    K = Kerr_Metric
    result = compute_curvature(K.spacetime())
    assert result.is_ricci_flat and result.ricci_scalar == 0
    # Gamma^r_tt = M Delta (r^2 - a^2 cos^2 theta) / rho^6, a standard result for Boyer-Lindquist coordinates
    known = K.M * K.delta * (K.r**2 - K.a**2 * sp.cos(K.theta) ** 2) / K.rho2**3
    assert K.fast_simplify(result.christoffel[1][0][0] - known) == 0


def test_kerr_kretschmann_matches_the_known_closed_form_numerically():
    """Symbolically this takes too long, so evaluate the Riemann tensor at sample points instead."""
    import itertools
    K = Kerr_Metric
    gamma = compute_curvature(K.spacetime()).christoffel
    n, X = 4, K.coords
    d_gamma = {(r, nu, s, m): sp.diff(gamma[r][nu][s], X[m])
               for r, nu, s, m in itertools.product(range(n), repeat=4)}

    def kretschmann_at(point):
        num = lambda expr: sp.N(expr.subs(point), 30)
        gam = [[[num(gamma[l][i][j]) for j in range(n)] for i in range(n)] for l in range(n)]
        g, g_inv = K.g.subs(point).evalf(30), K.g_inv.subs(point).evalf(30)
        mixed = {}
        for r, s, m, nu in itertools.product(range(n), repeat=4):
            if m < nu:                                    # R^r_{s m nu}, antisymmetric in (m, nu)
                mixed[(r, s, m, nu)] = (num(d_gamma[(r, nu, s, m)]) - num(d_gamma[(r, m, s, nu)])
                                        + sum(gam[r][m][l] * gam[l][nu][s] - gam[r][nu][l] * gam[l][m][s]
                                              for l in range(n)))

        def lowered(a, s, m, nu):
            if m == nu:
                return 0
            sign, (m, nu) = (1, (m, nu)) if m < nu else (-1, (nu, m))
            return sign * sum(g[a, r] * mixed[(r, s, m, nu)] for r in range(n))

        comps = {idx: lowered(*idx) for idx in itertools.product(range(n), repeat=4)}
        comps = {idx: v for idx, v in comps.items() if abs(v) > 1e-25}
        return sum(g_inv[a, p] * g_inv[b, q] * g_inv[c, r] * g_inv[d, s] * v * w
                   for (a, b, c, d), v in comps.items() for (p, q, r, s), w in comps.items())

    for mass, spin, radius, cos_theta in [(1, sp.Rational(1, 2), 4, sp.Rational(3, 5)),
                                          (2, 1, 7, sp.Rational(-1, 3))]:
        point = {K.M: mass, K.a: spin, K.r: radius, K.theta: sp.acos(cos_theta)}
        rho2 = radius**2 + spin**2 * cos_theta**2
        known = (48 * mass**2 * (radius**2 - spin**2 * cos_theta**2)
                 * (rho2**2 - 16 * radius**2 * spin**2 * cos_theta**2) / rho2**6)
        assert abs(kretschmann_at(point) - sp.N(known, 30)) < 1e-20


# --- the program itself ------------------------------------------------------------

def run(argv, typed=()):
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err), mock.patch("builtins.input", side_effect=list(typed)):
        try:
            code = project.main(argv)
        except SystemExit as exit_:
            code = exit_.code
    return code, out.getvalue(), err.getvalue()


def test_importing_the_modules_prints_nothing():
    """The original scripts did all their work (and printing) the moment they were imported."""
    out = io.StringIO()
    with redirect_stdout(out):
        import importlib
        for module in (Schwarzschild_Metric, FLRW_Metric, Kerr_Metric, project):
            importlib.reload(module)
    assert out.getvalue() == ""


def test_cli_schwarzschild():
    code, out, _ = run(["schwarzschild", "--kretschmann"])
    assert code == 0 and "Ricci-flat" in out and "48*M**2/r**6" in out


def test_cli_flrw_with_options():
    code, out, _ = run(["flrw", "--k", "0", "--scale-factor", "exp(H*t)", "--brief"])
    assert code == 0 and "12*H**2" in out and "Gamma^" not in out


def test_cli_latex_output():
    code, out, _ = run(["schwarzschild", "--latex", "--brief"])
    assert code == 0 and "\\frac" in out


def test_cli_menu():
    code, out, _ = run([], typed=["x", "9", "1"])
    assert code == 0 and "Please type one of the numbers" in out and "Schwarzschild metric" in out


def test_cli_quit_from_menu():
    code, out, _ = run([], typed=["q"])
    assert code == 1 and "Cancelled" in out


def test_cli_bad_arguments_give_clear_errors():
    assert run(["nonsense"])[0] == 2
    assert run(["schwarzschild", "--k", "0"])[0] == 2
    assert run(["flrw", "--k", "7"])[0] == 2
    code, _, err = run(["flrw", "--scale-factor", "exp((("])
    assert code == 2 and "couldn't read the scale factor" in err


def test_kerr_kretschmann_request_is_explained_not_ignored():
    code, out, _ = run(["kerr", "--kretschmann", "--brief"])
    assert code == 0 and "not computed for Kerr" in out


if __name__ == "__main__":
    tests = [(name, func) for name, func in sorted(globals().items()) if name.startswith("test_")]
    failed = 0
    for name, func in tests:
        try:
            func()
            print(f"ok    {name}")
        except Exception as error:                                   # noqa: BLE001
            failed += 1
            print(f"FAIL  {name}: {type(error).__name__}: {error}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)
