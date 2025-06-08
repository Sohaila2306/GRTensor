GENERAL RELATIVITY TENSOR CALCULUS IN PYTHON
=============================================

This repository provides symbolic computations of spacetime curvature for major metrics in general relativity using Python's SymPy library. It serves as a free and open-source alternative to Maple's GRTensor II.

PROJECT STRUCTURE
------------------
GeneralRelativity_Project.py      -> Master script for computing and displaying tensor quantities
Schwarzschild_Metric.py          -> Schwarzschild metric tensor and tensor operations
Kerr_Metric.py                   -> Kerr metric tensor (rotating spacetime)
FLRW_Metric.py                   -> FLRW metric tensor (cosmological spacetime)
README.txt                       -> Project documentation

WHAT THIS PROJECT DOES
-----------------------
- Defines major spacetime metrics:
  • Schwarzschild: spherical symmetry
  • Kerr: axial symmetry (rotating black hole)
  • FLRW: homogeneous and isotropic (cosmology)

- Computes:
  • Metric tensor (g_mu_nu)
  • Christoffel symbols (Gamma^lambda_mu_nu)
  • Ricci tensor (R_mu_nu)
  • Ricci scalar (R)

- Supports extensions: Riemann tensor, Kretschmann scalar, and geodesic equations.

------
Run the main script to compute tensors:
    python GeneralRelativity_Project.py

Output includes metrics and their Christoffel symbols, Ricci tensor, and scalar curvature. Kerr metric skips heavy computation by default.

NOTES
------
- All computations are symbolic and ideal for theoretical physics.
- Computationally expensive for rotating spacetimes like Kerr.
- Easily extensible for other metrics or additional tensor quantities.

BACKGROUND
-----------
Inspired by the functionality of GRTensor II (http://grtensor.phy.queensu.ca/), originally built for Maple. This project replicates similar capabilities using Python and SymPy.

