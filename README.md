# Topological Dark Solitons and the Core–Cusp Problem

This repository contains the theory, code, and data supporting:

> **Topological Dark Solitons from Mutual Chern–Simons Locking:  
> A Stable Resolution of the Core–Cusp Problem**

by **Xicplice**.

---

## Overview

We study a dark-sector model with two ultra-weakly coupled U(1) gauge fields linked by a mutual Chern–Simons term.  
On a topologically wound background this interaction generates a scalar mode $\varepsilon(\mathbf{x})\$ obeying a screened Helmholtz equation:

```math
\left(\nabla^2 - \frac{1}{\lambda_\varepsilon^2}\right)\varepsilon
= -\frac{8\pi G}{c^2}\rho_b(\mathbf{x})
```

with screening length:

```math
\lambda_\varepsilon = \frac{R}{\sqrt{n_1 n_2}}.
```

The exact spherical solution:

```math
\varepsilon(r) = \varepsilon_\infty\left(1 - e^{-r/\lambda_\varepsilon}\right)
```

induces a soliton-like baryonic density profile:

```math
\rho_b(r) = \frac{c^2\varepsilon_\infty}{8\pi G\lambda_\varepsilon^2}
\left(1 - \frac{2\lambda_\varepsilon}{r}e^{-r/\lambda_\varepsilon}\right),
```

with total defect mass:

```math
\Delta M = -1\,\frac{c^2\varepsilon_\infty\lambda_\varepsilon}{G}.
```

For $\lambda_\varepsilon \sim 10\ \mathrm{kpc}$, this removes
$\sim 10^{11} M_\odot$ from the inner halo—naturally converting an NFW cusp into a core while preserving large-scale ΛCDM phenomenology.

---

<table>
<tr>
<td><img src="plots/epsilon_profile.png" width="350"></td>
<td><img src="plots/density_profiles.png" width="350"></td>
</tr>
<tr>
<td><img src="plots/mass_profile.png" width="350"></td>
<td><img src="plots/rotation_curve.png" width="350"></td>
</tr>
</table>

---

## Repository Contents

- `paper.tex` — main LaTeX source of the paper  
- `docs/theory.md` — theoretical derivations and model explanation  
- `docs/numerics.md` — numerical methods and implementation details  
- `src/helmholtz_solver.py` — analytic Helmholtz–soliton solver with regularisation  
- `src/density_profiles.py` — soliton density, mass, and velocity tools  
- `src/locking_energy.py` — topological phase-locking (Josephson-like) sector  
- `src/stability_tests.py` — regression and consistency tests  
- `examples/milky_way_profile.py` — Milky Way–like soliton demonstration  
- `examples/parameter_scan.py` — parameter scan across $\lambda_\varepsilon\$ and $\varepsilon_\infty\$

---

## Installation

Clone and create a virtual environment:

```bash
git clone https://github.com/xicplice/topological-dark-solitons.git
cd topological-dark-solitons

python -m venv venv
source venv/bin/activate      # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

(Optional) Install the package in editable mode:

```bash
pip install -e .
```

---

# Quick Start

## Example 1: Milky Way–like soliton profile

From the repository root:

```bash
python -m examples.milky_way_profile
```

This prints:

- soliton parameters  
- example scalar field values  
- regularised density values  
- approximate defect mass $\Delta M\$

---

## Example 2: Parameter grid scanning

Explore how the defect mass scales with $\lambda_\varepsilon\$ and $\varepsilon_\infty\$:

```bash
python -m examples.parameter_scan \
    --R 1.0 \
    --n1 12 --n2 93 \
    --eps-inf-min 1e-7 --eps-inf-max 1e-5 --n-eps 5 \
    --lambda-factor-min 0.5 --lambda-factor-max 2.0 --n-lambda 5
```

This outputs a table:

```
lambda_factor   lambda_eps   eps_inf   defect_mass
```

---

# Development & Testing

A basic numerical sanity-check suite is included.

Run it using:

```bash
python -m src.stability_tests
```

This performs:

- finiteness and monotonicity checks of $\rho(r)\$  
- negative-defect-mass verification  
- rotation-curve smoothness checks  

Typical output:

```
=== Stability test suite ===
[OK] Test passed.
[OK] Test passed.
[OK] Test passed.

Overall status: PASS
```

---

# Citation

A `CITATION.cff` file is provided.  
GitHub natively supports this; citation metadata is exported automatically.

If you use this work, please cite it using the provided citation file.

---

# License

- **Code:** Apache-2.0  
- **Paper:** CC-BY-4.0  

See repository metadata for details.
