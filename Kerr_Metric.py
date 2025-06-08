import sympy as sp
from sympy import symbols, sin, cos, sqrt, simplify, Matrix

# Step 1: Define coordinates and parameters
t, r, theta, phi = symbols('t r theta phi')
M = symbols('M', real=True, positive=True)  # Mass parameter
a = symbols('a', real=True)  # Rotation parameter (angular momentum per mass)

# Step 2: Define auxiliary functions for Kerr metric
rho2 = r**2 + (a * cos(theta))**2
delta = r**2 - 2*M*r + a**2

# Step 3: Define Kerr metric components in Boyer-Lindquist coordinates
# The metric is a 4x4 symmetric matrix g_{mu nu}
g = Matrix(4, 4, lambda i, j: 0)

# g_tt
g[0, 0] = -(1 - (2*M*r)/rho2)

# g_tphi and g_phit (off-diagonal terms)
g[0, 3] = - (2*M*r*a*sin(theta)**2) / rho2
g[3, 0] = g[0, 3]

# g_rr
g[1, 1] = rho2 / delta

# g_thetatheta
g[2, 2] = rho2

# g_phiphi
g[3, 3] = (r**2 + a**2 + (2*M*r*a**2*sin(theta)**2)/rho2) * sin(theta)**2

# Coordinates list
coords = [t, r, theta, phi]

# Step 4: Print the Kerr metric neatly
def neat_matrix_print(matrix):
    for row in matrix.tolist():
        print("[", end=" ")
        for elem in row:
            elem_str = str(simplify(elem))
            print(f"{elem_str:>30}", end=" ")
        print("]")

print("\nKerr Metric tensor g (Boyer-Lindquist coordinates):")
neat_matrix_print(g)

# Optional: To extend, reuse Christoffel and Ricci tensor functions from Project I, but beware
# the calculations are much more complex and computationally heavy for Kerr.

print("\nNote: Christoffel symbols and Ricci tensor calculations for Kerr metric can be done similarly but require significant compute time.")
