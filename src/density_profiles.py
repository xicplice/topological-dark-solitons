"""
Convenience layer around :mod:`helmholtz_solver`.

Provides named imports for the most-used functions and a helper
for cumulative mass profiles.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from .helmholtz_solver import (
    SolitonParameters,
    epsilon_analytic,
    rho_analytic,
    rho_regularised,
    defect_mass,
    rotation_curve,
)


__all__ = [
    "SolitonParameters",
    "epsilon_analytic",
    "rho_analytic",
    "rho_regularised",
    "defect_mass",
    "rotation_curve",
    "mass_profile",
]


def mass_profile(
    r: np.ndarray,
    p: SolitonParameters,
    rho_halo: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> np.ndarray:
    """
    Cumulative mass profile M(r) for the soliton plus optional halo.

    This is similar to :func:`rotation_curve` but returns M(r)
    rather than v_c(r).

    Parameters
    ----------
    r : array_like
        Radii at which to evaluate the cumulative mass.
    p : SolitonParameters
        Soliton configuration.
    rho_halo : callable, optional
        Additional density component ρ_halo(r), if any.

    Returns
    -------
    np.ndarray
        Cumulative mass M(r).
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
    dr = r_int[1] - r_int[0]
    M_cum = np.cumsum(integrand) * dr

    return np.interp(r, r_int, M_cum)
