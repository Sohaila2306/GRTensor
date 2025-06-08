# GeneralRelativity_Project.py

import sympy as sp
from Schwarzschild_Metric import t as t_s, r as r_s, theta as theta_s, phi as phi_s, M as M_s, g as g_s, \
    compute_christoffel, compute_ricci
from Kerr_Metric import t as t_k, r as r_k, theta as theta_k, phi as phi_k, M as M_k, a as a_k, g as g_k, \
    neat_matrix_print as neat_print_kerr
from FLRW_Metric import t as t_f, r as r_f, theta as theta_f, phi as phi_f, a as a_f, k as k_f, g as g_f, \
    compute_christoffel as compute_christoffel_flrw, compute_ricci as compute_ricci_flrw, \
    neat_matrix_print as neat_print_flrw


def neat_matrix_print(matrix):
    """Generic neat matrix print, fallback."""
    for row in matrix.tolist():
        print("[", end=" ")
        for elem in row:
            elem_str = str(sp.simplify(elem))
            print(f"{elem_str:>30}", end=" ")
        print("]")


def compute_and_display(metric_name, g, coords):
    print(f"\n==== {metric_name} ====")
    # Compute Christoffel symbols
    Gamma = compute_christoffel(g, coords)
    print("\nMetric tensor g:")
    neat_matrix_print(g)

    print("\nChristoffel Symbols Γ^l_{ij} (nonzero terms):")
    n = len(coords)
    for l in range(n):
        for i in range(n):
            for j in range(n):
                if Gamma[l][i][j] != 0:
                    print(f"Γ^{l}_{i}{j} =", Gamma[l][i][j])

    # Compute Ricci tensor
    R = compute_ricci(Gamma, coords)
    print("\nRicci Tensor R_mu_nu:")
    neat_matrix_print(R)

    # Compute Ricci scalar
    g_inv = g.inv()
    R_scalar = sum(g_inv[i, j] * R[i, j] for i in range(n) for j in range(n))
    R_scalar = sp.simplify(R_scalar)
    print("\nRicci Scalar R:")
    print(R_scalar)


def main():
    # Schwarzschild Metric
    coords_s = [t_s, r_s, theta_s, phi_s]
    print("\n*** Schwarzschild Metric ***")
    compute_and_display("Schwarzschild Metric", g_s, coords_s)

    # Kerr Metric (Print only metric as Christoffel/Ricci heavy)
    print("\n*** Kerr Metric ***")
    print("Metric tensor g:")
    neat_print_kerr(g_k)
    print(
        "\nNote: Christoffel symbols and Ricci tensor calculations for Kerr metric can be done similarly but are computationally intensive.")

    # FLRW Metric
    coords_f = [t_f, r_f, theta_f, phi_f]
    print("\n*** FLRW Metric ***")
    compute_and_display("FLRW Metric", g_f, coords_f)


if __name__ == "__main__":
    main()
