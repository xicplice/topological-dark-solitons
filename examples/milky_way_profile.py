"""
Example: Milky Way-like soliton profile.

Chooses parameters such that the screening length λ_ε is of order
R / sqrt(n1 n2), and prints defect mass and sample values.
"""

from __future__ import annotations

import numpy as np

from src.helmholtz_solver import (
    SolitonParameters,
    epsilon_analytic,
    rho_regularised,
    defect_mass,
)


def main() -> None:
    R = 1.0  # code units
    n1, n2 = 12, 93
    eps_inf = 1e-6

    params = SolitonParameters(R=R, n1=n1, n2=n2, eps_inf=eps_inf, r_core=0.05)

    print("=== Soliton parameters ===")
    print(f"n1 = {n1}, n2 = {n2}, n1*n2 = {n1 * n2}")
    print(f"lambda_eps = {params.lambda_eps:.4e} (in units of R)")
    print(f"eps_inf    = {params.eps_inf:.4e}")
    print(f"r_core     = {params.r_core:.4e}")

    r = np.linspace(0.0, 5.0 * params.lambda_eps, 200)
    eps_r = epsilon_analytic(r, params)
    rho_r = rho_regularised(r, params)

    print("\nSample radii and fields:")
    for rr, ee, rh in zip(r[::40], eps_r[::40], rho_r[::40]):
        print(f"r = {rr:.3e}, eps = {ee:.3e}, rho = {rh:.3e}")

    dM = defect_mass(params)
    print(f"\nApproximate defect mass ΔM ≈ {dM:.4e} (code units).")


if __name__ == "__main__":
    main()
