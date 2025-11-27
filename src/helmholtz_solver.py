"""
Radial Helmholtz solver and soliton utilities.

Implements the scalar sector

    (∇² - λ^{-2}) ε(r) = - (8πG / c²) ρ_b(r)

with both the exact analytic profile and a finite-core regularisation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

# Physical constants (SI); can be set to 1 for code units.
G_SI = 6.67430e-11  # m^3 kg^-1 s^-2
C_SI = 299_792_458.0  # m s^-1


@dataclass
class SolitonParameters:
    """
    Container for soliton and screening parameters.

    Parameters
    ----------
    R : float
        Macroscopic scale (same units as lambda_eps).
    n1, n2 : int
        Topological integers.
    eps_inf : float
        Asymptotic value ε∞.
    lambda_eps : float or None
        Screening length λ_ε. If None, λ_ε = R / sqrt(n1 n2).
    r_core : float
        Core radius r_0 for regularisation.
    use_SI : bool
        If True: use physical G, c. If False: use G = c = 1
        (code units).
    """

    R: float
    n1: int
    n2: int
    eps_inf: float
    lambda_eps: Optional[float] = None
    r_core: float = 0.0
    use_SI: bool = False

    def __post_init__(self) -> None:
        if self.lambda_eps is None:
            if self.n1 <= 0 or self.n2 <= 0:
                raise ValueError(
                    "n1 and n2 must be positive to infer lambda_eps."
                )
            self.lambda_eps = self.R / np.sqrt(self.n1 * self.n2)

        if self.lambda_eps <= 0:
            raise ValueError("lambda_eps must be positive.")
        if self.r_core < 0:
            raise ValueError("core radius r_core must be non-negative.")

    @property
    def G(self) -> float:
        """Gravitational coupling (either SI or code units)."""
        return G_SI if self.use_SI else 1.0

    @property
    def c(self) -> float:
        """Speed of light (either SI or code units)."""
        return C_SI if self.use_SI else 1.0

    @property
    def rho_infinity(self) -> float:
        """
        Asymptotic density ρ∞ implied by the analytic profile.

        ρ∞ = c⁴ ε∞ / (8π G λ²).
        """
        lam = self.lambda_eps
        return (self.c ** 2 * self.eps_inf) / (8.0 * np.pi * self.G * lam**2)


def epsilon_analytic(r: np.ndarray, p: SolitonParameters) -> np.ndarray:
    """Analytic ε(r) = ε∞ (1 - exp(-r/λ)) profile."""
    r = np.asarray(r, dtype=float)
    lam = p.lambda_eps
    return p.eps_inf * (1.0 - np.exp(-r / lam))


def rho_analytic(r: np.ndarray, p: SolitonParameters) -> np.ndarray:
    r"""
    Exact density profile implied by the analytic ε(r):

        ρ_b(r) = (c² ε∞ / (8π G λ²)) [ 1 - (2λ / r) exp(-r/λ) ].

    The expression inherits a 1/r singularity at r → 0. For numerical
    purposes we avoid a literal division by zero and instead treat r=0
    as a placeholder value; in most applications you will want to use
    :func:`rho_regularised` instead.
    """
    r = np.asarray(r, dtype=float)
    lam = p.lambda_eps
    pref = (p.c ** 2 * p.eps_inf) / (8.0 * np.pi * p.G * lam**2)

    # Avoid exact division by zero; the singularity is regularised
    # separately.
    safe_r = np.where(r == 0.0, 1.0, r)
    return pref * (1.0 - (2.0 * lam / safe_r) * np.exp(-r / lam))


def rho_regularised(r: np.ndarray, p: SolitonParameters) -> np.ndarray:
    """
    Regularised density with a finite constant-density core.

    For r < r_core: ρ = ρ_core = ρ_analytic(r_core);
    For r ≥ r_core: ρ = ρ_analytic(r).

    If r_core == 0, this is identical to :func:`rho_analytic`.
    """
    r = np.asarray(r, dtype=float)
    rho = rho_analytic(r, p)

    if p.r_core > 0.0:
        rho_core = float(rho_analytic(np.array([p.r_core]), p)[0])
        mask = r < p.r_core
        rho[mask] = rho_core

    return rho


def defect_mass(
    p: SolitonParameters,
    r_max_factor: float = 20.0,
    n_points: int = 10_000,
) -> float:
    """
    Compute the total defect mass ΔM by integrating (ρ - ρ∞).

    ΔM = ∫ [ρ(r) - ρ∞] 4π r² dr,
    where ρ is the regularised density profile and ρ∞ is the
    asymptotic density at large r.

    Parameters
    ----------
    p : SolitonParameters
        Soliton configuration.
    r_max_factor : float
        Integration radius in units of λ_ε (default: 20).
    n_points : int
        Number of radial grid points.

    Returns
    -------
    float
        Approximate defect mass ΔM.
    """
    lam = p.lambda_eps
    r_max = r_max_factor * lam
    r = np.linspace(0.0, r_max, n_points)

    rho_reg = rho_regularised(r, p)
    rho_inf = p.rho_infinity

    integrand = (rho_reg - rho_inf) * 4.0 * np.pi * r**2
    return float(np.trapz(integrand, r))


def rotation_curve(
    r: np.ndarray,
    p: SolitonParameters,
    rho_halo: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> np.ndarray:
    """
    Compute circular velocity v_c(r) from soliton + optional halo:

        M(r) = ∫⁰ʳ 4π r'² [ρ_soliton + ρ_halo] dr' ,
        v_c²(r) = G M(r) / r.

    Parameters
    ----------
    r : array_like
        Radii at which to evaluate the circular velocity.
    p : SolitonParameters
        Soliton configuration.
    rho_halo : callable, optional
        Function rho_halo(r) giving any additional density component.

    Returns
    -------
    np.ndarray
        Circular velocities v_c(r).
    """
    r = np.asarray(r, dtype=float)
    if r.ndim == 0:
        r = r[None]

    r_int = np.linspace(0.0, float(r.max()), max(2_000, r.size * 5))
    rho_sol = rho_regularised(r_int, p)

    if rho_halo is not None:
        rho_tot = rho_sol + rho_halo(r_int)
    else:
        rho_tot = rho_sol

    integrand = 4.0 * np.pi * r_int**2 * rho_tot
    # Simple cumulative integral (left Riemann sum)
    dr = r_int[1] - r_int[0]
    M_cum = np.cumsum(integrand) * dr

    # Interpolate M(r) onto the requested radii.
    M_of_r = np.interp(r, r_int, M_cum)

    with np.errstate(divide="ignore", invalid="ignore"):
        v2 = p.G * M_of_r / np.where(r == 0.0, np.nan, r)
    v2[~np.isfinite(v2)] = 0.0

    return np.sqrt(v2)
