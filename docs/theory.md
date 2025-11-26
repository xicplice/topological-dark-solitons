# Theory: Topological Dark Solitons and Screened Helmholtz Profiles

This document sketches the theoretical background behind the code in
`src/`.

## 1. Two-sector gauge–scalar model

We consider two ultra-weakly coupled $U(1)$ gauge fields $A_\mu$ and
$B_\mu$ in the dark sector, linked by a mutual Chern–Simons interaction
on a topologically non-trivial background. Schematically,
\[
\mathcal{L} \supset -\frac{1}{4}F_{\mu\nu}F^{\mu\nu}
                   -\frac{1}{4}G_{\mu\nu}G^{\mu\nu}
                   + \kappa\,\epsilon^{\mu\nu\rho\sigma}
                     A_\mu \partial_\nu B_{\rho\sigma}
                   + \cdots.
\]

On a background with winding numbers $(n_1,n_2)$, the effective
low-energy dynamics can be recast in terms of a scalar mode
$\varepsilon(\mathbf{x})$ that captures the locked configuration of
the two sectors.

## 2. Screened Helmholtz equation

For non-relativistic sources with baryonic density $\rho_b(\mathbf{x})$
one arrives at the screened Helmholtz equation
\[
\biggl(\nabla^2 - \frac{1}{\lambda_\varepsilon^2}\biggr)\varepsilon
= -\frac{8\pi G}{c^4} \rho_b,
\]
where the screening length is set by the macroscopic scale $R$ and the
topological charges $n_1,n_2$,
\[
\lambda_\varepsilon = \frac{R}{\sqrt{n_1 n_2}}.
\]

In spherical symmetry the equation becomes
\[
\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d\varepsilon}{dr}\right)
 - \frac{\varepsilon}{\lambda_\varepsilon^2}
 = -\frac{8\pi G}{c^4}\rho_b(r).
\]

## 3. Analytic soliton profile

For a broad class of sources the scalar profile can be written
\[
\varepsilon(r) = \varepsilon_\infty
\left(1 - e^{-r/\lambda_\varepsilon}\right),
\]
which solves the homogeneous Helmholtz equation at large $r$ and
approaches $\varepsilon_\infty$ asymptotically.

The corresponding source density implied by the equation of motion is
\[
\rho_b(r) = \frac{c^4\varepsilon_\infty}{8\pi G\lambda_\varepsilon^2}
\left(1 - \frac{2\lambda_\varepsilon}{r}e^{-r/\lambda_\varepsilon}\right).
\]
This profile has a $1/r$ singularity at the origin but rapidly
approaches a constant at $r \gtrsim \lambda_\varepsilon$.

The total defect mass,
\[
\Delta M = \int d^3\mathbf{x}\,\bigl[\rho_b(r) - \rho_\infty\bigr],
\]
is negative for the parameter ranges of interest and can be shown to
scale as
\[
\Delta M = -2\frac{c^4\varepsilon_\infty\lambda_\varepsilon}{G}.
\]

## 4. Core regularisation

To render the profile finite at the origin we replace the central
region $r < r_0$ by a constant density equal to $\rho_b(r_0)$. This
is implemented in `rho_regularised` in `src/helmholtz_solver.py`:

- for $r < r_0$, $\rho(r) = \rho_b(r_0)$,
- for $r \ge r_0$, $\rho(r) = \rho_b(r)$.

This regularisation preserves the exterior solution but removes the
singularity.

## 5. Topological phase locking

The same topological integers $(n_1,n_2)$ control a Josephson-like
locking energy between the two sectors,
\[
E(\alpha) = E_0\bigl[1 - \cos(n_1 n_2 \alpha)\bigr],
\]
where $\alpha$ is a relative phase. For small angles,
\[
E(\alpha) \simeq \frac{1}{2}K_{\rm eff}\,\alpha^2,
\qquad
K_{\rm eff} = E_0(n_1n_2)^2,
\]
which defines a stiffness scaling quadratically with $n_1n_2$.

This sector is implemented in `src/locking_energy.py` and is useful
for estimating phase-locking energy scales and small-angle dynamics.

## 6. Connection to the core–cusp problem

For $\lambda_\varepsilon$ of order $10\ \mathrm{kpc}$ the negative
defect mass in the interior of a halo is of order
$10^{11}M_\odot$, sufficient to transform an NFW cusp into a core
in Milky Way-sized galaxies while leaving large-scale structure and
CMB constraints essentially unchanged.

The code in this repository provides ready-to-use tools to compute:

- density profiles (regularised and unregularised),
- cumulative masses and circular velocities,
- basic sanity checks on the numerical behaviour of the profiles.
