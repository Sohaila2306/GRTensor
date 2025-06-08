import sympy as sp
from sympy import symbols, diag, sin, simplify, Matrix

# Step 1: Define the coordinates of spacetime and the mass parameter
# We use 't' for time, 'r' for radial distance, and 'theta' and 'phi' for angles.
t, r, theta, phi = symbols('t r theta phi')
M = symbols('M', real=True, positive=True)  # Mass of the central object (like a black hole)

# Step 2: Define the Schwarzschild metric tensor, which describes how spacetime is curved
# around a spherical, non-rotating mass.
# This is a diagonal matrix with four components (time, radial, and two angular parts).
g = diag(
    -(1 - 2 * M / r),  # Time-time component (how time is affected by gravity)
    1 / (1 - 2 * M / r),  # Radial component (how radial distances behave)
    r ** 2,  # Theta-theta component (related to the polar angle)
    r ** 2 * sin(theta) ** 2  # Phi-phi component (related to the azimuthal angle)
)

coords = [t, r, theta, phi]  # Keep track of the order of coordinates


# Step 3: Calculate the Christoffel symbols, which capture how vectors change
# when moved around in curved spacetime. They’re essential for describing geodesics
# and curvature.
def compute_christoffel(metric, coordinates):
    n = len(coordinates)
    metric_inv = metric.inv()  # Inverse of the metric tensor
    Gamma = [[[0] * n for _ in range(n)] for _ in range(n)]  # Initialize Christoffel symbols as zeros

    for l in range(n):
        for i in range(n):
            for j in range(n):
                # Sum over k to compute the Christoffel symbol component Γ^l_ij
                val = 0
                for k in range(n):
                    val += metric_inv[l, k] * (
                            sp.diff(metric[k, j], coordinates[i]) +
                            sp.diff(metric[k, i], coordinates[j]) -
                            sp.diff(metric[i, j], coordinates[k])
                    )
                Gamma[l][i][j] = simplify(0.5 * val)  # Symmetrize and simplify the result
    return Gamma


Gamma = compute_christoffel(g, coords)


# Step 4: Compute the Ricci tensor, which tells us about the curvature caused by matter
# and energy in spacetime. It’s derived from the Christoffel symbols.
def compute_ricci(Gamma, coordinates):
    n = len(coordinates)
    R = Matrix.zeros(n, n)  # Initialize Ricci tensor with zeros

    for mu in range(n):
        for nu in range(n):
            total = 0
            for alpha in range(n):
                # Calculate derivatives and sums over the Christoffel symbols to find R_mu_nu
                total += (sp.diff(Gamma[alpha][mu][nu], coordinates[alpha]) -
                          sp.diff(Gamma[alpha][mu][alpha], coordinates[nu]) +
                          sum(Gamma[alpha][mu][beta] * Gamma[beta][nu][alpha] -
                              Gamma[alpha][nu][beta] * Gamma[beta][mu][alpha]
                              for beta in range(n)))
            R[mu, nu] = simplify(total)
    return R


R = compute_ricci(Gamma, coords)

# Step 5: Calculate the Ricci scalar by contracting the Ricci tensor with the inverse metric
# It’s a single number representing overall curvature of spacetime.
g_inv = g.inv()
R_scalar = sum(g_inv[i, j] * R[i, j] for i in range(4) for j in range(4))
R_scalar = simplify(R_scalar)


# Helper function: Nicely print matrices row by row for better readability
def neat_matrix_print(matrix):
    for row in matrix.tolist():
        print("[", end=" ")
        for elem in row:
            # Convert each element to a simplified string and align it for neatness
            elem_str = str(sp.simplify(elem))
            print(f"{elem_str:>20}", end=" ")
        print("]")


# Final outputs: show the metric tensor, Christoffel symbols, Ricci tensor, and Ricci scalar

print("\nMetric tensor g:")
neat_matrix_print(g)

print("\nChristoffel Symbols Γ^l_{ij}:")
for l in range(4):
    for i in range(4):
        for j in range(4):
            if Gamma[l][i][j] != 0:
                print(f"Γ^{l}_{i}{j} =", Gamma[l][i][j])

print("\nRicci Tensor R_mu_nu:")
neat_matrix_print(R)

print("\nRicci Scalar R:")
print(R_scalar)

print(
    "\nNote: The Kretschmann scalar, a measure of overall spacetime curvature, can be computed using more advanced tools like EinsteinPy or by fully expanding the Riemann tensor.")
