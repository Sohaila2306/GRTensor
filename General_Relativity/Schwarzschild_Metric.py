"""Schwarzschild_Metric.py - the spacetime outside a spherical, non-rotating mass.

Found by Karl Schwarzschild in 1916, a few weeks after Einstein published the
field equations. It describes the gravity of a star, a planet or a (non-rotating)
black hole, in coordinates where r is the circumference of a sphere divided by 2*pi.

    ds^2 = -(1 - 2M/r) dt^2 + dr^2 / (1 - 2M/r) + r^2 (dtheta^2 + sin^2(theta) dphi^2)

Units are G = c = 1, so M is a length: the event horizon sits at r = 2M.

Run this file directly to see the curvature:    python Schwarzschild_Metric.py
"""

import sympy as sp

from Curvature_Tools import Spacetime, compute_curvature, print_report

t, r, theta, phi = sp.symbols("t r theta phi", real=True)
M = sp.symbols("M", positive=True)          # mass of the central object

coords = [t, r, theta, phi]

# f(r) shows up in both the time and radial parts. It vanishes at the horizon r = 2M.
f = 1 - 2 * M / r
g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta) ** 2)


def spacetime() -> Spacetime:
    return Spacetime(
        name="Schwarzschild metric",
        metric=g,
        coords=coords,
        description="Vacuum outside a spherical, non-rotating mass M (G = c = 1).",
    )


if __name__ == "__main__":
    print_report(compute_curvature(spacetime(), with_kretschmann=True))
