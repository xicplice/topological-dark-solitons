# Numerics: Helmholtz Solver and Stability Checks

This document describes the numerical approach implemented in the
`src/` package.

## 1. Design philosophy

The goal is to keep the numerics simple, transparent, and easily
portable into other codes:

- All computations use `numpy` arrays.
- No external ODE/PDE solvers are required for the current analytic
  profile; numerical integrals are performed with simple quadrature
  (`np.trapz`) and cumulative sums.
- Units are controlled by a flag (`use_SI`) inside `SolitonParameters`:
  - in code units (`G = c = 1`),
  - or in physical SI units.

## 2. Scalar sector implementation

The main scalar and density utilities live in
`src/helmholtz_solver.py`:

- `SolitonParameters`: dataclass containing $(R, n_1, n_2,
  \varepsilon_\infty, \lambda_\varepsilon, r_{\rm core})$ and unit
  choices.
- `epsilon_analytic(r, p)`: analytic scalar profile
  $\varepsilon(r)$.
- `rho_analytic(r, p)`: implied density profile $\rho_b(r)$ from the
  Helmholtz equation.
- `rho_regularised(r, p)`: finite-core regularisation.
- `defect_mass(p, ...)`: numerical integral of the excess density
  over a finite radial domain.
- `rotation_curve(r, p, rho_halo=None)`: circular velocity from the
  soliton plus an optional additional halo component.

Computational details:

- The radial integrals use evenly spaced grids in $r$:
  - For `defect_mass`, the grid covers $[0, r_{\max}]$ with
    $r_{\max} \sim 20\,\lambda_\varepsilon$ by default.
  - For `rotation_curve`, an internal grid is built up to
    `r.max()`, with sufficient resolution to ensure smooth
    interpolation.

- `rho_regularised` uses a vectorised mask:
  - compute the unregularised profile everywhere,
  - identify indices with $r < r_{\rm core}$,
  - replace them by a constant value obtained from
    $\rho_{\rm analytic}(r_{\rm core})$.

## 3. Density and mass profiles

`src/density_profiles.py` provides a convenience layer:

- Re-export of the main scalar and density functions.
- Helper functions for cumulative mass `mass_profile` and circular
  velocity wrappers useful in fitting and parameter scans.

These utilities avoid duplication of logic and keep the physics
“single-sourced” in `helmholtz_solver.py`.

## 4. Stability and regression tests

`src/stability_tests.py` implements basic numerical consistency checks:

- Detect NaNs and infinities in `rho_regularised` on a wide radial
  grid.
- Verify that the defect mass has the expected (negative) sign for
  typical parameters.
- Check smoothness of the rotation curve via finite differences.

These tests are intended for development-time validation rather than
a full physical stability analysis.

To run the tests:

```bash
python -m src.stability_tests
