"""GeneralRelativity_Project.py - curvature of the Schwarzschild, Kerr and FLRW spacetimes.

For each spacetime this works out, symbolically, the Christoffel symbols, the
Ricci tensor and the Ricci scalar (and the Kretschmann scalar where that is
practical).

    python GeneralRelativity_Project.py                  # pick from a menu
    python GeneralRelativity_Project.py schwarzschild
    python GeneralRelativity_Project.py kerr --brief     # skip the long Christoffel list
    python GeneralRelativity_Project.py flrw --k 0 --scale-factor "exp(H*t)"
    python GeneralRelativity_Project.py all --kretschmann
"""

from __future__ import annotations

import argparse
import sys

import sympy as sp

import FLRW_Metric
import Kerr_Metric
import Schwarzschild_Metric
from Curvature_Tools import compute_curvature, print_report

CHOICES = {
    "schwarzschild": "Schwarzschild - spherical, non-rotating mass (a black hole or a star)",
    "kerr": "Kerr - rotating black hole (takes ~10 seconds)",
    "flrw": "FLRW - an expanding, homogeneous universe",
}

KERR_KRETSCHMANN_NOTE = (
    "Note: the Kretschmann scalar is not computed for Kerr (it takes far too long symbolically).\n"
    "      Its known closed form is in the README."
)


def build_spacetime(name: str, k_value=None, scale_factor=None):
    if name == "schwarzschild":
        return Schwarzschild_Metric.spacetime()
    if name == "kerr":
        return Kerr_Metric.spacetime()
    return FLRW_Metric.spacetime(k_value, scale_factor)


def run(name: str, kretschmann: bool = False, latex: bool = False, brief: bool = False,
        k_value=None, scale_factor=None) -> None:
    spacetime = build_spacetime(name, k_value, scale_factor)
    want_k = kretschmann and name != "kerr"
    print(f"\nComputing {spacetime.name}...", flush=True)
    result = compute_curvature(
        spacetime, with_kretschmann=want_k,
        progress=lambda stage: print(f"  - {stage}", flush=True))
    print_report(result, latex=latex, show_christoffel=not brief)
    if kretschmann and name == "kerr":
        print("\n" + KERR_KRETSCHMANN_NOTE)


def choose_from_menu() -> str:
    names = list(CHOICES)
    print("Which spacetime would you like to explore?\n")
    for number, name in enumerate(names, start=1):
        print(f"  {number} - {CHOICES[name]}")
    print(f"  {len(names) + 1} - All three")
    while True:
        answer = input(f"\nEnter 1-{len(names) + 1} (or q to quit): ").strip().lower()
        if answer in ("q", "quit", "exit"):
            raise KeyboardInterrupt
        if answer.isdigit() and 1 <= int(answer) <= len(names):
            return names[int(answer) - 1]
        if answer == str(len(names) + 1) or answer == "all":
            return "all"
        print("  Please type one of the numbers above.")


def build_parser() -> argparse.ArgumentParser:
    def curvature_sign(text):
        if text not in ("-1", "0", "1"):
            raise argparse.ArgumentTypeError("must be -1, 0 or 1")
        return int(text)

    parser = argparse.ArgumentParser(
        description="Symbolic curvature (Christoffel symbols, Ricci tensor, Ricci scalar) "
                    "for the Schwarzschild, Kerr and FLRW spacetimes. "
                    "Run with no arguments for a menu.")
    parser.add_argument("metric", nargs="?", choices=[*CHOICES, "all"],
                        help="which spacetime to analyse")
    parser.add_argument("--kretschmann", action="store_true",
                        help="also compute the Kretschmann scalar (not available for Kerr)")
    parser.add_argument("--brief", action="store_true",
                        help="don't list the Christoffel symbols")
    parser.add_argument("--latex", action="store_true",
                        help="print results as LaTeX instead of plain text")
    flrw = parser.add_argument_group("FLRW options")
    flrw.add_argument("--k", dest="k_value", type=curvature_sign, metavar="{-1,0,1}",
                      help="fix the spatial curvature (default: leave k symbolic)")
    flrw.add_argument("--scale-factor", metavar="EXPR",
                      help='fix the scale factor as a function of t, e.g. "exp(H*t)" or "t**(2/3)" '
                           "(default: leave a(t) as an unknown function)")
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    scale_factor = None
    if args.scale_factor is not None:
        try:
            scale_factor = sp.sympify(args.scale_factor, locals={"t": FLRW_Metric.t})
        except (sp.SympifyError, SyntaxError, TypeError):
            parser.error(f"couldn't read the scale factor {args.scale_factor!r}")
    if (args.k_value is not None or scale_factor is not None) and args.metric not in (None, "flrw"):
        parser.error("--k and --scale-factor only apply to the flrw metric")

    try:
        name = args.metric or choose_from_menu()
        for item in (list(CHOICES) if name == "all" else [name]):
            run(item, args.kretschmann, args.latex, args.brief, args.k_value, scale_factor)
    except (KeyboardInterrupt, EOFError):
        print("\nCancelled.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
