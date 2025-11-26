# Topological Dark Solitons and the Core–Cusp Problem

This repository contains the theory, code, and data supporting the paper

> **Topological Dark Solitons from Mutual Chern–Simons Locking:  
> A Stable Resolution of the Core–Cusp Problem**

by **Xicplice**.

## Overview

We study a dark-sector model with two ultra-weakly coupled $U(1)$
gauge fields linked by a mutual Chern–Simons term. On a topologically
wound background this interaction generates a scalar mode
$\varepsilon(\mathbf{x})$ obeying a screened Helmholtz equation
\[
\left(\nabla^2 - \frac{1}{\lambda_\varepsilon^2}\right)\varepsilon
= -\frac{8\pi G}{c^4}\rho_b(\mathbf{x}),
\]
with screening length
\[
\lambda_\varepsilon = \frac{R}{\sqrt{n_1 n_2}}.
\]

The exact spherical solution
\[
\varepsilon(r) = \varepsilon_\infty\bigl(1 - e^{-r/\lambda_\varepsilon}\bigr)
\]
corresponds to a negative-mass soliton with density
\[
\rho_b(r) = \frac{c^4\varepsilon_\infty}{8\pi G\lambda_\varepsilon^2}
\left(1 - \frac{2\lambda_\varepsilon}{r}e^{-r/\lambda_\varepsilon}\right),
\]
and total defect mass
\[
\Delta M = -2\frac{c^4\varepsilon_\infty\lambda_\varepsilon}{G}.
\]

For $\lambda_\varepsilon \sim 10\ \mathrm{kpc}$ this removes
$\sim 10^{11}M_\odot$ from the inner halo, converting an
NFW cusp into a core without modifying large-scale $\Lambda$CDM.

## Repository contents

- `paper.tex`: main LaTeX source of the paper.
- `docs/theory.md`: pedagogical explanation of the model and derivations.
- `docs/numerics.md`: numerical methods and convergence tests.
- `src/helmholtz_solver.py`: radial Helmholtz solver with core regularisation.
- `src/density_profiles.py`: analytic and numerical density / mass profiles.
- `src/locking_energy.py`: phase-locking energy, stiffness, and scaling.
- `src/stability_tests.py`: perturbation experiments and diagnostics.
- `examples/milky_way_profile.py`: example Milky Way-like profile and diagnostics.
- `examples/parameter_scan.py`: parameter study over $(\lambda_\varepsilon, \varepsilon_\infty)$.

## Installation

Clone and create a virtual environment:

```bash
git clone https://github.com/YOUR_USERNAME/topological-dark-solitons.git
cd topological-dark-solitons

python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

pip install -r requirements.txt
