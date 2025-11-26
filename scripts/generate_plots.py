#!/usr/bin/env python3
"""
Generate scientific plots for README and Zenodo using the
model equations implemented in the repository.

This script uses:
- epsilon_analytic(r)
- rho_analytic(r)
- rho_regularised(r)
- defect_mass integration
- rotation_curve(r)

from src/helmholtz_solver.py

All figures are saved into: plots/
"""

import numpy as np
import matplotlib.pyplot as plt
import os

from src.helmholtz_solver import (
    SolitonParameters,
    epsilon_analytic,
    rho_analytic,
    rho_regularised,
    defect_mass,
    rotation_curve,
)

# ================================
# 1. Create output directory
# ================================
OUTDIR = "plots"
os.makedirs(OUTDIR, exist_ok=True)

# ================================
# 2. Define physical model parameters
# ================================
R = 1.0
n1, n2 = 12, 93
eps_inf = 1e-6
r_core = 0.05

params = SolitonParameters(
    R=R,
    n1=n1,
    n2=n2,
    eps_inf=eps_inf,
    r_core=r_core,
    use_SI=False  # code units for all plots
)

lam = params.lambda_eps

print("=== Plot Generation Parameters ===")
print(f"R = {R}")
print(f"n1, n2 = {n1}, {n2}")
print(f"n1*n2 = {n1*n2}")
print(f"lambda_eps = {lam}")
print(f"eps_inf = {eps_inf}")
print(f"r_core = {r_core}")
print("===================================")

# Radial grid
r = np.linspace(0, 5 * lam, 2000)
r_safe = np.where(r == 0, 1e-12, r)

# ================================
# 3. Compute model functions
# ================================
eps_r = epsilon_analytic(r, params)
rho_a = rho_analytic(r_safe, params)
rho_r = rho_regularised(r_safe, params)

# Cumulative mass used in rotation curve:
# (rotation_curve internally does the integral on finer grid)
# But here we compute M(r) explicitly for plotting:
integrand = 4 * np.pi * r_safe**2 * rho_r
mass_r = np.cumsum(integrand) * (r_safe[1] - r_safe[0])

v_c = rotation_curve(r, params)

# ================================
# 4. Plot 1: Scalar field profile ε(r)
# ================================
plt.figure(figsize=(6, 4))
plt.plot(r, eps_r, linewidth=2)
plt.xlabel("r")
plt.ylabel("ε(r)")
plt.title("Analytic Scalar Field Profile ε(r)")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/epsilon_profile.png", dpi=300)
plt.close()

# ================================
# 5. Plot 2: Density Profiles
# ================================
plt.figure(figsize=(6, 4))
plt.plot(r, rho_a, label="Analytic ρ(r)", linewidth=2)
plt.plot(r, rho_r, label="Regularised ρ(r)", linewidth=2)
plt.xlabel("r")
plt.ylabel("ρ(r)")
plt.title("Density Profiles (Analytic vs Regularised)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/density_profiles.png", dpi=300)
plt.close()

# ================================
# 6. Plot 3: Cumulative Mass M(r)
# ================================
plt.figure(figsize=(6, 4))
plt.plot(r, mass_r, linewidth=2)
plt.xlabel("r")
plt.ylabel("M(r)")
plt.title("Cumulative Defect Mass Profile M(r)")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/mass_profile.png", dpi=300)
plt.close()

# ================================
# 7. Plot 4: Rotation Curve
# ================================
plt.figure(figsize=(6, 4))
plt.plot(r, v_c, linewidth=2)
plt.xlabel("r")
plt.ylabel("v_c(r)")
plt.title("Rotation Curve from Soliton Model")
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rotation_curve.png", dpi=300)
plt.close()

# ================================
# 8. Print completion message
# ================================
print("\nPlots generated successfully in: plots/")
print("Files:")
print(" - epsilon_profile.png")
print(" - density_profiles.png")
print(" - mass_profile.png")
print(" - rotation_curve.png")
