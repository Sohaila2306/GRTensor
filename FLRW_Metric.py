import sympy as sp
from sympy import symbols, diag, sin, simplify, Matrix

# -------------------------
# Step 1: Define coordinates and parameters
# -------------------------
t, r, theta, phi = symbols('t r theta phi', real=True)
a = sp.Function('a')(t)  # Scale factor a(t)
k = symbols('k', real=True)  # Spatial curvature constant: k = -1,0,1

# -------------------------
# Step 2: Define the FLRW metric
# -------------------------
# Spatial part depends on k
# Define the radial function depending on curvature k
def f_r(r, k):
    if k == 0:
        return r
    elif k == 1:
        return sp.sin(r)
    elif k == -1:
        return sp.sinh(r)
    else:
        # General formula often written as 1/(1 - k r^2)
        return r

# For symbolic, define the metric elements as general functions of k and r:
g_tt = -1
g_rr = a**2 / (1 - k * r**2)
g_thth = a**2 * r**2
g_phph = a**2 * r**2 * sin(theta)**2

g = diag(g_tt, g_rr, g_thth, g_phph)
coords = [t, r, theta, phi]

# -------------------------
# Step 3: Compute Christoffel symbols
# -------------------------
def compute_christoffel(g, coords):
    n = len(coords)
    g_inv = g.inv()
    Gamma = [[[0]*n for j in range(n)] for i in range(n)]
    for l in range(n):
        for i in range(n):
            for j in range(n):
                sum_term = 0
                for k_ in range(n):
                    sum_term += g_inv[l, k_] * (
                        sp.diff(g[k_, j], coords[i]) +
                        sp.diff(g[k_, i], coords[j]) -
                        sp.diff(g[i, j], coords[k_])
                    )
                Gamma[l][i][j] = simplify(1/2 * sum_term)
    return Gamma

Gamma = compute_christoffel(g, coords)

# -------------------------
# Step 4: Compute Ricci tensor
# -------------------------
def compute_ricci(Gamma, coords):
    n = len(coords)
    R = Matrix.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            sum_term = 0
            for alpha in range(n):
                sum_term += (
                    sp.diff(Gamma[alpha][mu][nu], coords[alpha]) -
                    sp.diff(Gamma[alpha][mu][alpha], coords[nu]) +
                    sum(
                        Gamma[alpha][mu][beta] * Gamma[beta][nu][alpha] -
                        Gamma[alpha][nu][beta] * Gamma[beta][mu][alpha]
                        for beta in range(n)
                    )
                )
            R[mu, nu] = simplify(sum_term)
    return R

R = compute_ricci(Gamma, coords)

# -------------------------
# Step 5: Compute Ricci scalar
# -------------------------
g_inv = g.inv()
R_scalar = sum([g_inv[i, j] * R[i, j] for i in range(4) for j in range(4)])
R_scalar = simplify(R_scalar)

# -------------------------
# Helper function to print matrices neatly
# -------------------------
def neat_matrix_print(M):
    for row in M.tolist():
        print("[", end=" ")
        for elem in row:
            elem_str = str(sp.simplify(elem))
            print(f"{elem_str:>20}", end=" ")
        print("]")

# -------------------------
# Final output
# -------------------------
print("\nFLRW Metric tensor g:")
neat_matrix_print(g)

print("\nChristoffel Symbols Γ^l_{ij} (nonzero terms):")
for l in range(4):
    for i in range(4):
        for j in range(4):
            if Gamma[l][i][j] != 0:
                print(f"Γ^{l}_{i}{j} =", Gamma[l][i][j])

print("\nRicci Tensor R_mu_nu:")
neat_matrix_print(R)

print("\nRicci Scalar R:")
print(R_scalar)

print("\nNote: This symbolic output can be specialized by substituting scale factor a(t) and curvature k.")
