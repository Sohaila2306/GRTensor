# General Relativity Tensor Calculus in Python

Symbolic curvature calculations for three of the most important spacetimes in general relativity: **Schwarzschild** (a non-rotating black hole), **Kerr** (a rotating one) and **FLRW** (an expanding universe). For each one the code works out the metric's Christoffel symbols, Ricci tensor and Ricci scalar exactly, as formulas, using [SymPy](https://www.sympy.org). Where it is practical it also computes the Kretschmann scalar.

It is a free, open-source alternative for the kind of work people used to do in Maple with [GRTensor II](http://grtensor.phy.queensu.ca/). You can use it to check homework, explore a metric you've invented, or just watch how Einstein's equations behave.

```
$ python GeneralRelativity_Project.py schwarzschild --kretschmann

Christoffel symbols (non-zero, Gamma^l_{ij} = Gamma^l_{ji} so only i <= j is listed):
  Gamma^t_(t,r) = M/(r*(-2*M + r))
  Gamma^r_(t,t) = M*(-2*M + r)/r**3
  Gamma^r_(r,r) = M/(r*(2*M - r))
  ...
Ricci tensor R_{mu nu} (non-zero components):
  all components are zero (Ricci-flat: a vacuum solution)

Ricci scalar R = 0
Kretschmann scalar K = 48*M**2/r**6
```

## Why curvature? A short history

General relativity says that gravity is not a force pulling things through space but the geometry of spacetime itself, and that matter and energy decide how that geometry curves. Making that idea precise took mathematics that was developed decades before anyone knew what to use it for.

- **1854.** Bernhard Riemann gives his famous lecture on the foundations of geometry, introducing the idea that space can be curved in ways you can measure from the inside, without stepping outside it.
- **1869.** Elwin Christoffel works out the symbols that bear his name. They tell you how the coordinate directions twist and stretch from point to point, which is what you need to compare vectors in a curved space.
- **Around 1900.** Gregorio Ricci-Curbastro and his student Tullio Levi-Civita develop the "absolute differential calculus", today's tensor calculus, and the Ricci tensor is named after the former.
- **1912-1915.** Einstein has the physical idea (the equivalence principle) but struggles with the maths. His friend Marcel Grossmann points him to Riemann, Christoffel and Ricci. In November 1915 Einstein presents the final field equations and shows they explain the unexplained 43 arcseconds per century in the orbit of Mercury.
- **1916.** Karl Schwarzschild, serving in the German army on the Russian front, finds the first exact solution of the equations within weeks of their publication. He dies a few months later.
- **1917-1935.** Einstein applies the equations to the whole universe. Alexander Friedmann (1922) and Georges Lemaitre (1927) find solutions describing an expanding universe. Edwin Hubble's observations (1929) show the expansion is real, and Howard Robertson and Arthur Walker (1935) prove that the metric must take this form for any universe that looks the same everywhere.
- **1963.** After almost half a century of searching, Roy Kerr, a New Zealand mathematician, finds the solution for a rotating mass. Robert Boyer and Richard Lindquist introduce the coordinates used today in 1967.
- **1960s onward.** Hand calculation of these tensors is long and easy to get wrong, so relativists were among the first people to use computer algebra. Systems like GRTensor (written for Maple at Queen's University in Canada), xAct for Mathematica, Cadabra, and SageManifolds followed.

### What the tensors mean

| Quantity | What it tells you |
|---|---|
| **Metric** `g` | How to measure distances and time intervals in these coordinates |
| **Christoffel symbols** `Gamma` | How the coordinate grid itself changes from point to point. They appear in the equation for the path of a freely falling particle, which is where "gravity" shows up |
| **Riemann tensor** | The full curvature. It describes tidal forces: how two nearby falling particles drift together or apart |
| **Ricci tensor** `R_mu_nu` | The part of the curvature that responds to matter and energy. It describes how a small ball of freely falling dust changes volume. In Einstein's equations it is tied directly to the matter present |
| **Ricci scalar** `R` | The Ricci tensor boiled down to one number at each point |
| **Kretschmann scalar** `K` | The "size" of the full curvature as one number that doesn't depend on your choice of coordinates, so it can tell a real singularity from a coordinate artefact |

In empty space Einstein's equations reduce to `R_mu_nu = 0`. Both Schwarzschild and Kerr are exactly such solutions, which makes "is the Ricci tensor zero?" one of the best correctness checks there is for code like this.

## The three spacetimes

### Schwarzschild: a non-rotating black hole (or star)

The geometry outside any spherical, non-rotating mass `M`. Everything simple about gravity in the Solar System is a consequence of it: the bending of starlight seen in the 1919 eclipse, the precession of Mercury's orbit, the radar delay of signals passing the Sun, and the clock corrections that GPS satellites need (about 38 microseconds a day). Taken to its limit, with the mass squeezed inside `r = 2M`, it describes a black hole.

It also teaches a good lesson. The metric looks singular at `r = 2M`, but the Kretschmann scalar there is perfectly finite, which shows that the horizon is an artefact of the coordinates (this took physicists decades to fully accept, until coordinates like Kruskal-Szekeres in 1960 made it plain). At `r = 0` the Kretschmann scalar `48M^2/r^6` blows up, and that singularity is real.

### Kerr: a rotating black hole

Real stars spin, and so do the black holes they collapse into. Kerr is the solution for a spinning mass with parameters `M` and `a = J/M`. Setting `a = 0` gives Schwarzschild back. It predicts effects that Schwarzschild can't: spacetime itself is dragged around the hole, there is an "ergosphere" where nothing can stay still, and, as Roger Penrose showed in 1969, energy can be extracted from the hole's rotation.

Astrophysical black holes are believed to be described by Kerr, which is why it matters so much today. It underlies the analysis of the gravitational waves first detected by LIGO in 2015 and the images of the shadows of M87* (2019) and of our own galaxy's Sagittarius A* (2022) from the Event Horizon Telescope.

### FLRW: an expanding universe

Friedmann, Lemaitre, Robertson and Walker's metric describes a universe that is the same in every place and every direction but changes with time. It has just two ingredients: the scale factor `a(t)`, which says how big the universe is, and the spatial curvature `k` (+1 closed, 0 flat, -1 open). Plug the Ricci tensor it produces into Einstein's equations and you get the Friedmann equations, the starting point of modern cosmology.

This is the geometry behind the Big Bang model, the cosmic microwave background found in 1965 by Penzias and Wilson, and the discovery in 1998 that the expansion is accelerating, which led to the 2011 Nobel Prize in Physics.

## What to expect

Every result here can be checked against a published value, and the tests do exactly that.

| Spacetime | Ricci tensor | Ricci scalar | Kretschmann scalar |
|---|---|---|---|
| Schwarzschild | 0 | 0 | `48 M^2 / r^6` |
| Kerr | 0 | 0 | `48 M^2 (r^2 - a^2 cos^2(theta)) (rho^4 - 16 r^2 a^2 cos^2(theta)) / rho^12`, with `rho^2 = r^2 + a^2 cos^2(theta)` (not computed symbolically; see below) |
| FLRW | `R_tt = -3 a''/a`, plus space-space components (see the program output) | `6 (a''/a + a'^2/a^2 + k/a^2)` | `12 [ (a''/a)^2 + ((a'^2 + k)/a^2)^2 ]` |

Special cases that make good sanity checks, all reproduced by the code:

- **de Sitter space** (`k = 0`, `a = exp(H t)`): `R = 12 H^2`, `K = 24 H^4`, and `R_mu_nu = 3 H^2 g_mu_nu`.
- **Matter-dominated universe** (`k = 0`, `a = t^(2/3)`): `R = 4 / (3 t^2)`, which equals `8 pi rho` for dust.

## Install

You need Python 3.9 or newer and SymPy:

```bash
pip install -r requirements.txt
```

## Use

Run the main script with no arguments and pick from a menu:

```bash
python GeneralRelativity_Project.py
```

Or go straight to a spacetime:

```bash
python GeneralRelativity_Project.py schwarzschild
python GeneralRelativity_Project.py kerr --brief
python GeneralRelativity_Project.py flrw
python GeneralRelativity_Project.py all
```

| Option | What it does |
|---|---|
| `--kretschmann` | Also compute the Kretschmann scalar (Schwarzschild and FLRW) |
| `--brief` | Skip the Christoffel symbols, which are long for Kerr |
| `--latex` | Print results as LaTeX, ready to paste into a paper |
| `--k {-1,0,1}` | FLRW only: fix the spatial curvature instead of leaving it as a symbol |
| `--scale-factor EXPR` | FLRW only: fix `a(t)`, for example `"exp(H*t)"` or `"t**(2/3)"` |

For example, to see that a matter-dominated universe has `R = 4/(3 t^2)`:

```bash
python GeneralRelativity_Project.py flrw --k 0 --scale-factor "t**(2/3)" --brief
```

Each metric file can also be run on its own (`python Kerr_Metric.py`), or imported:

```python
import sympy as sp
import FLRW_Metric
from Curvature_Tools import compute_curvature

H = sp.symbols("H", positive=True)
result = compute_curvature(FLRW_Metric.spacetime(k_value=0, scale_factor=sp.exp(H * FLRW_Metric.t)))
print(result.ricci_scalar)       # 12*H**2
```

### Try your own metric

Anything that is a metric matrix plus a list of coordinates will work. Here is a charged black hole (Reissner-Nordstrom):

```python
import sympy as sp
from Curvature_Tools import Spacetime, compute_curvature, print_report

t, r, theta, phi = sp.symbols("t r theta phi", real=True)
M, Q = sp.symbols("M Q", positive=True)

f = 1 - 2*M/r + Q**2/r**2
g = sp.diag(-f, 1/f, r**2, r**2 * sp.sin(theta)**2)

print_report(compute_curvature(Spacetime("Reissner-Nordstrom", g, [t, r, theta, phi]),
                               with_kretschmann=True))
```

The Ricci scalar comes out as zero (electromagnetic fields are traceless) while `R_tt` doesn't, because the electric field carries energy. The code gives `K = 8(6M^2 r^2 - 12 M Q^2 r + 7 Q^4)/r^8`, matching the published result. The functions don't assume four dimensions, so a 2-sphere or a 5-dimensional metric works too.

## How long does it take?

Timings measured on the machine used to test this project (yours will differ a little):

| Spacetime | Time |
|---|---|
| Schwarzschild (including Kretschmann) | about 1 second |
| FLRW (including Kretschmann) | about 1 second |
| Kerr (Christoffel symbols, Ricci tensor and scalar) | about 8 seconds |

Kerr used to be considered too heavy to compute this way, because SymPy's general `simplify` can take many minutes on it. The trick in `Kerr_Metric.py` is to avoid asking SymPy to simplify trigonometry. It writes the inverse metric in closed form and treats `sin(theta)` and `cos(theta)` as ordinary variables `S` and `C` with `S^2 + C^2 = 1`, then uses fast polynomial arithmetic. The Kretschmann scalar for Kerr is still too slow with this approach, so it is not computed in the program. Its closed form in the table above was checked numerically at sample points in the tests instead.

## Checking the results

```bash
python test_general_relativity.py     # no extra packages needed
# or, if you use pytest:
python -m pytest
```

The tests (about 20 seconds) check, among other things:

- Schwarzschild's Christoffel symbols against their textbook forms, and that its Ricci tensor is exactly zero
- Kerr's Ricci tensor is exactly zero, its inverse metric really is the inverse, it reduces to Schwarzschild when `a = 0`, and `Gamma^r_tt` matches the known formula
- FLRW against the Friedmann-equation ingredients, de Sitter space and the dust universe
- A 2-sphere, whose Ricci tensor must be `g_ij / a^2`
- That importing any module prints nothing and runs no heavy calculation
- The command-line options, including bad input

## Project files

| File | What it is |
|---|---|
| `GeneralRelativity_Project.py` | The main program: menu, options and printing |
| `Curvature_Tools.py` | The tensor calculus: Christoffel symbols, Ricci, Riemann, Kretschmann |
| `Schwarzschild_Metric.py` | The Schwarzschild metric |
| `Kerr_Metric.py` | The Kerr metric, with its fast simplifier |
| `FLRW_Metric.py` | The FLRW metric, with optional `k` and scale factor |
| `test_general_relativity.py` | Tests |
| `requirements.txt` | Dependencies |

## Limitations

- Everything is symbolic and exact, which is wonderful for understanding and slow for large problems. For numbers (orbits, ray tracing) you want a numerical code.
- Units are `G = c = 1` and the signature is `(-, +, +, +)`.
- Geodesic equations and the Einstein tensor are not included yet. They are natural next steps, since the Christoffel symbols and Ricci tensor are already there.
- Results aren't automatically simplified into the prettiest possible form; two correct expressions can look different. Use `sp.simplify` on the difference to compare them.

## Further reading

- Sean Carroll, *Spacetime and Geometry*: a readable modern introduction
- Misner, Thorne and Wheeler, *Gravitation*: the classic reference, and the source of this project's sign conventions
- Robert Wald, *General Relativity*: a more mathematical treatment
- The GRTensor II documentation, to see what a mature tensor-algebra package offers
