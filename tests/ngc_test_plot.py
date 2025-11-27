import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from astropy import constants as const
from astropy import units as u
import os

# ================================================================
# 0. CONFIGURATION
# ================================================================

# List of galaxies. Add more entries here when you have more files.
# For each galaxy, we try the files in 'files' in order.
GALAXIES = [
    {
        "name": "NGC 2403",
        "files": ["NGC2403_rotmod.dat", "NGC2403.csv"],
        "R_fit_DM": 10.0,   # kpc, radius inside which we compare ΔM
    },
    # Example for future:
    # {
    #     "name": "NGC 3198",
    #     "files": ["NGC3198_rotmod.dat", "NGC3198.csv"],
    #     "R_fit_DM": 15.0,
    # },
]

# Cosmological prior on rho_vac (very loose, tweakable)
RHO_VAC_MIN = 1e-3  # Msun/pc^3
RHO_VAC_MAX = 1e-1  # Msun/pc^3
SIGMA_LOG_RHO_VAC = 0.5  # controls how strong the prior is

# Weight of mass-deficit penalty relative to rotation-curve chi^2
ALPHA_MASS_DEFICIT = 0.5   # tweak: 0 = ignore ΔM, >1 = enforce strongly


# ================================================================
# 1. DATA LOADING
# ================================================================

def load_sparc_like(path):
    """
    Load SPARC .dat or CSV in standard format:
    Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul
    """
    if path.endswith(".dat"):
        df = pd.read_csv(
            path,
            comment="#",
            sep=r"\s+",
            names=["Rad","Vobs","errV","Vgas","Vdisk","Vbul","SBdisk","SBbul"],
            engine="python"
        )
    else:
        df = pd.read_csv(path)

    # Ensure columns exist
    for col in ["Rad","Vobs","errV","Vgas","Vdisk","Vbul"]:
        if col not in df.columns:
            df[col] = 0.0

    # Convert to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Extract & clean
    r   = df["Rad"].values
    v   = df["Vobs"].values
    dv  = df["errV"].fillna(5.0).values
    vg  = df["Vgas"].values
    vd  = df["Vdisk"].values
    vbul= df["Vbul"].values

    mask = np.isfinite(r) & np.isfinite(v)
    return {
        "r":   r[mask],
        "v":   v[mask],
        "dv":  dv[mask],
        "vg":  vg[mask],
        "vd":  vd[mask],
        "vbul":vbul[mask],
    }

# Load all galaxies
gal_data = []
for g in GALAXIES:
    fname = None
    for cand in g["files"]:
        if os.path.exists(cand):
            fname = cand
            break
    if fname is None:
        raise FileNotFoundError(f"No file found for galaxy {g['name']}")
    print(f"Loading {g['name']} from {fname}")
    data = load_sparc_like(fname)
    data["name"] = g["name"]
    data["R_fit_DM"] = g["R_fit_DM"]
    gal_data.append(data)

N_gal = len(gal_data)
print(f"\nLoaded {N_gal} galaxies.")


# ================================================================
# 2. SOLITON + NFW MODEL
# ================================================================

G_KPC = 4.301e-6  # G in kpc (km/s)^2 Msun^-1
C_KMS = 2.998e5   # c in km/s

def dm_velocity_profile(r_points, eps, lam, rho0, rs):
    r_points = np.asarray(r_points)
    r_grid = np.logspace(-2, np.log10(r_points.max()*1.2), 800)

    rho_vac = (C_KMS**2 * eps) / (8*np.pi*G_KPC*lam**2)
    rho_nfw = rho0 / ((r_grid/rs) * (1 + r_grid/rs)**2)
    rho_sol = rho_vac * (1 - (2*lam/r_grid)*np.exp(-r_grid/lam))
    rho_dm  = rho_nfw + rho_sol - rho_vac
    rho_dm[rho_dm < 0] = 0.0

    M_dm = cumulative_trapezoid(4*np.pi*r_grid**2 * rho_dm, r_grid, initial=0)
    v_dm = np.sqrt(G_KPC * M_dm / r_grid)
    return np.interp(r_points, r_grid, v_dm)

def enclosed_DM_mass(R_kpc, eps, lam, rho0, rs):
    """Enclosed DM mass inside radius R_kpc (numeric)."""
    r_grid = np.logspace(-2, 3, 2000)  # up to 1e3 kpc
    rho_vac = (C_KMS**2 * eps) / (8*np.pi*G_KPC*lam**2)
    rho_nfw = rho0 / ((r_grid/rs)*(1 + r_grid/rs)**2)
    rho_sol = rho_vac * (1 - (2*lam/r_grid)*np.exp(-r_grid/lam))
    rho_dm  = rho_nfw + rho_sol - rho_vac
    rho_dm[rho_dm < 0] = 0.0
    M_dm = cumulative_trapezoid(4*np.pi*r_grid**2 * rho_dm, r_grid, initial=0)
    return np.interp(R_kpc, r_grid, M_dm)  # Msun


# ================================================================
# 3. GLOBAL PARAMETER VECTOR
# ================================================================
# PARAMS:
#   p[0] = lam (kpc, global)
#   p[1] = log10(eps) (global)
#   for each galaxy i:
#       p[2 + 3*i + 0] = log10(rho0_i)
#       p[2 + 3*i + 1] = log10(rs_i)
#       p[2 + 3*i + 2] = log10(ML_disk_i)  (M/L relative factor; 0 => 1.0)

def unpack_params(p):
    lam = p[0]
    eps = 10**p[1]
    per_gal = []
    for i in range(N_gal):
        base = 2 + 3*i
        log_rho0 = p[base + 0]
        log_rs   = p[base + 1]
        log_ML   = p[base + 2]
        per_gal.append({
            "rho0": 10**log_rho0,
            "rs":   10**log_rs,
            "ML":   10**log_ML,
        })
    return lam, eps, per_gal


# ================================================================
# 4. CHI-SQUARED WITH PRIORS
# ================================================================

def total_chi2(p):
    lam, eps, per_gal = unpack_params(p)

    if lam <= 0 or eps <= 0:
        return 1e30

    # 4.1 rotation-curve chi^2
    chi2_rc = 0.0
    for g, pars in zip(gal_data, per_gal):
        r   = g["r"]
        v   = g["v"]
        dv  = g["dv"]
        vg  = g["vg"]
        vd0 = g["vd"]
        vbul= g["vbul"]

        rho0 = pars["rho0"]
        rs   = pars["rs"]
        ML   = pars["ML"]

        v_dm = dm_velocity_profile(r, eps, lam, rho0, rs)
        v_disk_scaled = vd0 * np.sqrt(ML)  # V ∝ sqrt(M)

        v_bary_sq = vg**2 + v_disk_scaled**2 + vbul**2
        v_tot = np.sqrt(v_bary_sq + v_dm**2)

        chi2_rc += np.sum(((v - v_tot)/dv)**2)

    # 4.2 cosmological prior on rho_vac
    lam_u = lam * u.kpc
    eps_inf = eps
    rho_vac = (const.c**2 * eps_inf / (8*np.pi*const.G * lam_u**2)).to(u.Msun/u.pc**3).value

    log_rho_vac = np.log10(rho_vac)
    log_mid = np.log10(np.sqrt(RHO_VAC_MIN*RHO_VAC_MAX))
    chi2_rho = ((log_rho_vac - log_mid)/SIGMA_LOG_RHO_VAC)**2

    # 4.3 mass-deficit penalties for each galaxy
    DeltaM_analytic = - (const.c**2 * eps_inf * lam_u / const.G).to(u.Msun).value

    chi2_DM = 0.0
    if DeltaM_analytic != 0.0:
        for g, pars in zip(gal_data, per_gal):
            Rfit = g["R_fit_DM"]
            rho0 = pars["rho0"]
            rs   = pars["rs"]

            M_with    = enclosed_DM_mass(Rfit, eps, lam, rho0, rs)
            M_without = enclosed_DM_mass(Rfit, 0.0, lam, rho0, rs)
            DeltaM_num = M_without - M_with  # Msun

            ratio = DeltaM_num / DeltaM_analytic
            # softly prefer ratio ~ 1, but don't blow up if sign flips
            chi2_DM += (ratio - 1.0)**2

    chi2_total = chi2_rc + chi2_rho + ALPHA_MASS_DEFICIT * chi2_DM
    return chi2_total


# ================================================================
# 5. RUN OPTIMISATION
# ================================================================

# Initial guesses
lam0 = 4.0
log_eps0 = np.log10(3e-7)

p0 = [lam0, log_eps0]
for g in gal_data:
    # crude per-galaxy guesses
    p0 += [
        np.log10(8e6),     # rho0
        np.log10(16.0),    # rs
        0.0,               # log10(ML_disk factor) = 1
    ]

# Bounds
bounds = []
bounds.append((0.5, 20.0))       # lam
bounds.append((-9.0, -4.0))     # log10 eps
for g in gal_data:
    bounds.append((5.0, 9.5))   # log rho0
    bounds.append((0.7, 2.0))   # log rs (5–100 kpc)
    bounds.append((-0.3, 0.3))  # log ML factor (~0.5–2)

print("\nRunning global optimisation...")
res = minimize(total_chi2, p0, bounds=bounds)
p_best = res.x
lam_best, eps_best, per_gal_best = unpack_params(p_best)

# compute final chi2 pieces to get reduced chi2 per galaxy
chi2_tot = total_chi2(p_best)
N_params = len(p_best)
N_points = sum(len(g["r"]) for g in gal_data)
chi2_red_global = chi2_tot / max(1, (N_points - N_params))


print("\n================ GLOBAL BEST-FIT ================")
print(f"lam (kpc) = {lam_best:.3f}")
print(f"eps       = {eps_best:.3e}")
print(f"rho_vac   = {(const.c**2*eps_best/(8*np.pi*const.G*(lam_best*u.kpc)**2)).to(u.Msun/u.pc**3):.3e}")
print(f"Global reduced chi^2 ≈ {chi2_red_global:.2f}")

DeltaM_analytic = - (const.c**2 * eps_best * (lam_best*u.kpc) / const.G).to(u.Msun)
print(f"Analytic ΔM = {DeltaM_analytic:.3e}")


# ================================================================
# 6. PLOTS + PER-GALAXY DIAGNOSTICS
# ================================================================

for g, pars in zip(gal_data, per_gal_best):
    r   = g["r"]
    v   = g["v"]
    dv  = g["dv"]
    vg  = g["vg"]
    vd0 = g["vd"]
    vbul= g["vbul"]

    rho0 = pars["rho0"]
    rs   = pars["rs"]
    ML   = pars["ML"]

    v_dm = dm_velocity_profile(r, eps_best, lam_best, rho0, rs)
    v_disk_scaled = vd0 * np.sqrt(ML)
    v_bary_sq = vg**2 + v_disk_scaled**2 + vbul**2
    v_tot = np.sqrt(v_bary_sq + v_dm**2)

    chi2_rc = np.sum(((v - v_tot)/dv)**2)
    dof = len(r) - 3  # per-galaxy free params: rho0, rs, ML (global lam,eps shared)
    chi2_red = chi2_rc / max(1, dof)

    # Plot
    plt.figure(figsize=(14,6))
    plt.errorbar(r, v, yerr=dv, fmt='ko', label="Observed Data")
    plt.plot(r, vg,  'c:', label="Gas")
    plt.plot(r, v_disk_scaled, 'y--', label="Stars (scaled)")
    plt.plot(r, v_dm,   'm--', label="Dark Matter (soliton+NFW)")
    plt.plot(r, v_tot,  'r-', linewidth=3,
             label=f"Total Fit (χ² ≈ {chi2_red:.2f})")
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Velocity (km/s)")
    plt.title(f"Full Rotation Curve Decomposition: {g['name']}")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    # Mass deficit diagnostic
    Rfit = g["R_fit_DM"]
    M_with    = enclosed_DM_mass(Rfit, eps_best, lam_best, rho0, rs)
    M_without = enclosed_DM_mass(Rfit, 0.0,     lam_best, rho0, rs)
    DeltaM_num = M_without - M_with

    ratio = DeltaM_num / DeltaM_analytic.value

    print(f"\n=== {g['name']} diagnostics ===")
    print(f"rho0     = {rho0:.3e} Msun/kpc^3")
    print(f"rs       = {rs:.3f} kpc")
    print(f"ML_disk  = {ML:.3f} (relative to template)")
    print(f"chi2_red (RC only) ≈ {chi2_red:.2f}")
    print(f"ΔM_numeric(<{Rfit} kpc) = {DeltaM_num:.3e} Msun")
    print(f"ΔM_numeric / ΔM_analytic = {ratio:.3f}")

    # simple quantisation proxy K^2 = (Rmax/lam)^2
    Rmax = r.max()
    K2 = (Rmax/lam_best)**2
    print(f"K^2 proxy = (Rmax/lam)^2 ≈ {K2:.2f}")
