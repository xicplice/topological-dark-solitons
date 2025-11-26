"""
Basic numerical stability and sanity tests for the soliton profiles.

This module is not a full physical stability analysis; instead it
provides regression-style checks that the implemented profiles are
finite, smooth, and have the expected sign for the defect mass.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .helmholtz_solver import (
    SolitonParameters,
    rho_regularised,
    defect_mass,
    rotation_curve,
)


@dataclass
class StabilityResult:
    ok: bool
    messages: list[str]


def _check_finiteness(p: SolitonParameters) -> StabilityResult:
    msgs: list[str] = []
    r = np.logspace(-4, 2, 2_000) * p.lambda_eps
    rho = rho_regularised(r, p)

    if not np.all(np.isfinite(rho)):
        msgs.append("Non-finite values detected in rho_regularised.")
    if np.any(rho < -1e6 * abs(p.rho_infinity)):
        msgs.append("Suspiciously large negative densities detected.")

    return StabilityResult(ok=not msgs, messages=msgs)


def _check_defect_mass_sign(p: SolitonParameters) -> StabilityResult:
    msgs: list[str] = []
    dm = defect_mass(p)

    # For the models of interest, we expect negative defect mass.
    if dm >= 0:
        msgs.append(f"Defect mass sign unexpected: ΔM = {dm:.3e} >= 0.")
    return StabilityResult(ok=not msgs, messages=msgs)


def _check_rotation_curve_smoothness(p: SolitonParameters) -> StabilityResult:
    msgs: list[str] = []
    r = np.linspace(1e-4 * p.lambda_eps, 20.0 * p.lambda_eps, 1_000)
    v = rotation_curve(r, p)

    if not np.all(np.isfinite(v)):
        msgs.append("Non-finite values detected in rotation_curve.")

    # Finite-difference check for wild oscillations.
    dv = np.diff(v)
    if np.any(np.abs(dv) > 10.0 * np.nanmax(np.abs(v))):
        msgs.append("Large jumps in rotation curve detected.")
    return StabilityResult(ok=not msgs, messages=msgs)


def run_all_tests(p: SolitonParameters) -> list[StabilityResult]:
    """
    Run all basic stability checks for a given parameter set.

    Returns a list of StabilityResult objects. A global pass occurs
    if all individual tests have result.ok == True.
    """
    tests = [
        _check_finiteness,
        _check_defect_mass_sign,
        _check_rotation_curve_smoothness,
    ]
    return [fn(p) for fn in tests]


def main() -> None:
    """
    Run the standard test suite for a representative parameter choice.
    """
    # Representative Milky Way-like configuration (code units).
    p = SolitonParameters(
        R=1.0,
        n1=12,
        n2=93,
        eps_inf=1e-6,
        r_core=0.05,
        use_SI=False,
    )

    results = run_all_tests(p)
    all_ok = all(res.ok for res in results)

    print("=== Stability test suite ===")
    for res in results:
        if res.ok:
            print("[OK] Test passed.")
        else:
            print("[FAIL]")
            for msg in res.messages:
                print("   -", msg)

    print("\nOverall status:", "PASS" if all_ok else "FAIL")


if __name__ == "__main__":
    main()
