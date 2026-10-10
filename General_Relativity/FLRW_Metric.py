"""FLRW_Metric.py - the Friedmann-Lemaitre-Robertson-Walker universe.

The metric of a universe that looks the same everywhere and in every direction
(homogeneous and isotropic) but can expand or shrink with time:

    ds^2 = -dt^2 + a(t)^2 [ dr^2 / (1 - k r^2) + r^2 (dtheta^2 + sin^2(theta) dphi^2) ]

a(t) is the scale factor and k is the curvature of space:
k = +1 closed (like a 3-sphere), k = 0 flat, k = -1 open. Units are c = 1.

By default a(t) is left as an unknown function and k as a symbol, so the answer
holds for every FLRW universe. Pass a number for k, and/or an explicit scale
factor, to look at a specific one:

    spacetime(k=0, scale_factor=sp.exp(H * t))        # de Sitter space

Run this file directly to see the curvature:    python FLRW_Metric.py
"""

import sympy as sp

from Curvature_Tools import Spacetime, compute_curvature, print_report

t, r, theta, phi = sp.symbols("t r theta phi", real=True)
k = sp.symbols("k", real=True)              # spatial curvature: -1, 0 or +1
a = sp.Function("a")(t)                     # scale factor a(t)

coords = [t, r, theta, phi]


def metric(k_value=None, scale_factor=None) -> sp.Matrix:
    """The FLRW metric, optionally with a fixed k and/or a chosen scale factor."""
    kk = k if k_value is None else sp.sympify(k_value)
    aa = a if scale_factor is None else sp.sympify(scale_factor)
    return sp.diag(
        -1,
        aa**2 / (1 - kk * r**2),
        aa**2 * r**2,
        aa**2 * r**2 * sp.sin(theta) ** 2,
    )


g = metric()


def spacetime(k_value=None, scale_factor=None) -> Spacetime:
    description = "Homogeneous, isotropic universe with scale factor a(t) and spatial curvature k (c = 1)."
    if k_value is not None or scale_factor is not None:
        description += f"  Here k = {k if k_value is None else k_value}, a(t) = {a if scale_factor is None else scale_factor}."
    return Spacetime(
        name="FLRW metric",
        metric=metric(k_value, scale_factor),
        coords=coords,
        description=description,
    )


if __name__ == "__main__":
    print_report(compute_curvature(spacetime(), with_kretschmann=True))
