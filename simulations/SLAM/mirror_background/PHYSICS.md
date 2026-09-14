# Physics guide: injected protons with plasma backgrounds

This package extends Lars's continuous proton source in a prescribed magnetic field. It is an early mirror-background study supporting the SLAM project; it is not a SLAM equilibrium or optimization result. The injected particles are already ionized protons. No neutral transport or ionization is simulated.

## What a background means in the code

A species is an ensemble of computational particles, each with a position, momentum, and weight. A weight states how many physical particles a computational particle represents. A thermal background has random velocities with zero mean drift; it is not motionless and does not behave like a prescribed fluid velocity field.

For a nonrelativistic Maxwellian, each velocity component has standard deviation

\[
v_{\mathrm{th},s}=\sqrt{\frac{e T_s[\mathrm{eV}]}{m_s}}.
\]

WarpX input momenta are normalized as \(u_j=\gamma v_j/c\). The Gaussian input spread is therefore `sqrt(T_s*q_e/m_s)/clight`. Density is in m^-3, position in metres, and time in seconds. `background_stopping.background_temperature` uses **kelvin**, so 10 eV must be converted to about 116,045 K for that parameter. Diagnostic momentum storage has its own openPMD units; analysis must read those units instead of assuming input normalization.

In the reservoir cases, `do_not_push = 1` holds the background ion positions fixed and disables their ordinary field pusher. Binary collisions can still change their momenta. This is a finite, spatially fixed ion bath, not a thermostat that restores a Maxwellian after every step. In the fully kinetic case, both electrons and background ions move under the particle pusher.

## Four model choices

| Variant | Background treatment | Beam interactions | Self-consistent electric field |
| --- | --- | --- | --- |
| `ion_background` | Spatially fixed thermal proton bath | Prescribed magnetic field; no collision operator | Disabled through charge deposition settings |
| `ion_coulomb` | Same ion bath | Magnetic orbit plus beam-ion binary Coulomb collisions | Disabled through charge deposition settings |
| `hybrid` | Same ion bath plus implicit Maxwellian electron reservoir | Beam-ion binary collisions plus analytic electron slowing | Disabled through charge deposition settings |
| Fully kinetic electrons and ions | Mobile thermal proton and electron populations | Binary Coulomb interactions among configured species pairs, electric forces, and prescribed magnetic field | Electrostatic Poisson solve with deposition enabled |

The first three cases use a selected bath density of 1e16 m^-3, ion temperature 1 eV, and, where applicable, electron temperature 10 eV. They deliberately omit collective charge-separation physics. Their density is a collision-model input and does not establish that a dense electron plasma is resolved on the mesh.

The fully kinetic demonstration uses physical electron/proton masses, `n_e = n_i = 1e10 m^-3`, `T_e = 10 eV`, and `T_i = 1 eV`. Both backgrounds use the same spatial mask and Cartesian `NUniformPerCell` positions with `2 2 2` particles per cell per species. Matching positions, density, and opposite charge cancel their initial deposited charge to rounding accuracy. Their thermal velocities differ and their later charge densities can separate.

The initialization window is `r < 0.75 m` and `0.25 < z < 4.75 m`. It keeps the short demonstration away from the longitudinal seam and most walls. It is a finite plasma column with edges, not a pre-established equilibrium. The successful runs use a grounded analytic bore of radius 0.9 m, approximating the STL inner surface; the faceted STL itself stalled the electrostatic solve. Check the generated EB mask and initial particle masks in the diagnostics; see RESULTS.md for the geometry comparison.

The continuously added proton beam contributes positive charge in the kinetic case. Equal initial electron and background-ion densities do not neutralize subsequently injected beam ions automatically. Lars's inherited flux of 1e5 m^-2 s^-1 and approximately 3e4 m/s normal speed imply a local beam density scale of only a few m^-3, enormously below this background density. Any beam-driven collective response will consequently be extremely small with that source strength.

## Fields and collisions are separate physics

For species \(s\), the orbit pusher implements the Lorentz force,

\[
\frac{d\mathbf x}{dt}=\mathbf v,\qquad
\frac{d\mathbf p}{dt}=q_s(\mathbf E+\mathbf v\times\mathbf B_{\mathrm{ext}}).
\]

The kinetic electrostatic model deposits \(\rho=\sum_s q_s n_s\), solves

\[
\nabla^2\phi=-\rho/\epsilon_0,\qquad \mathbf E=-\nabla\phi,
\]

and gathers that electric field back to the particles. This includes collective **electrostatic** effects. It does not evolve a self-consistent electromagnetic magnetic perturbation. The external HDF5 magnetic field is prescribed; its generating coils are not advanced by WarpX.

`pairwisecoulomb` is a stochastic binary-collision operator following the Perez et al. formulation. It groups and pairs computational particles locally and changes their momenta to model accumulated Coulomb scattering. It does not integrate the electromagnetic force from every microscopic encounter. Beam-ion scattering can change direction and exchange energy; collisions with electrons can also transfer beam energy. Intra-species pairs allow thermal populations to scatter internally.

The binary operators use the explicit comparison choice `CoulombLog = 10`. This is a model parameter, not a measured value or a calibration to SLAM. Unequal particle weights mean binary updates need not conserve total energy and momentum exactly without correction. The injected beam and bath weights are especially different here, so energy checks and correction warnings matter.

The hybrid electron operator instead applies the low-speed, Maxwellian-electron approximation

\[
\frac{d\mathbf v_b}{dt}=-\alpha\mathbf v_b,\qquad
\mathbf v_b(t+\Delta t)=\mathbf v_b(t)e^{-\alpha\Delta t},
\]

with \(\alpha\propto n_e\ln\Lambda/(m_b T_e^{3/2})\) for fixed beam charge and physical electron mass. It scales all momentum components together, so this operator alone reduces speed without rotating its direction. Orbital motion and ion scattering can separately change direction. It assumes \(m_b\gg m_e\) and beam speed much smaller than electron thermal speed; the inherited approximately 9-16 eV proton source with 10 eV electrons is in that low-speed regime.

`background_stopping` computes its own Coulomb logarithm from its bath parameters in the inspected source. Its logarithm is not overridden by the binary operators' `CoulombLog = 10`. It provides deterministic drag, not the complete electron velocity-diffusion operator or a kinetic electron population. Energy deposited into its implicit reservoir is not stored in electron particles.

**Do not add analytic electron drag to a kinetic beam-electron collision pair for the same background.** That would count electron slowing twice. The hybrid and kinetic descriptions are alternative models, not successive collision terms that all remain enabled.

## Why the kinetic density and timestep differ

The inherited 2 m by 2 m by 5 m box with a `32 32 80` grid has \(\Delta x=\Delta y=\Delta z=0.0625\) m. Relevant electron scales are

\[
\lambda_{De}=\sqrt{\frac{\epsilon_0 eT_e}{n_e e^2}},\quad
\omega_{pe}=\sqrt{\frac{n_e e^2}{\epsilon_0m_e}},\quad
\Omega_{ce}=\frac{e|B|}{m_e}.
\]

For the selected kinetic density, \(\lambda_{De}\approx0.235\) m and \(\omega_{pe}\approx5.64\times10^6\) rad/s. Including both initial thermal species gives

\[
\lambda_D^{-2}=\lambda_{De}^{-2}+\lambda_{Di}^{-2},
\qquad \lambda_D\approx0.0709\ \mathrm{m}.
\]

The combined length is only about 1.13 cells: it passes this basic scale check but is not a spatial-convergence demonstration. At 1e16 m^-3 with the same temperatures, these lengths are 1000 times smaller, so the present mesh cannot resolve that dense kinetic plasma with the same explicit setup. Simply entering the reservoir density into the kinetic input would not produce a trustworthy dense-plasma calculation.

The audited field maximum is approximately 0.00923865 T. It implies \(\Omega_{ce,\max}\approx1.625\times10^9\) rad/s. The kinetic smoke timestep is 5e-11 s, for which \(\Omega_{ce,\max}\Delta t\approx0.0812\). It also meets the selected conservative plasma-frequency and thermal-cell-crossing checks:

\[
\Delta t\leq\min\left(
\frac{0.1}{\Omega_{ce,\max}},
\frac{0.1}{\omega_{pe}},
\frac{0.2\Delta x}{4v_{\mathrm{th},e}}
\right).
\]

These are timestep design choices for this test, not universal accuracy guarantees. The finite Maxwellian sample can contain velocities beyond four thermal standard deviations; measured velocity extrema and fields must still be checked.

The reservoir tests use 1e-7 s steps for an initial 20 microsecond horizon. With physical protons this gives \(\Omega_{ci,\max}\Delta t\approx0.0885\). The kinetic smoke uses 200 steps for **10 nanoseconds = 0.01 microseconds**, with a half-timestep companion at the same physical horizon. That spans about 2.6 fastest electron gyroperiods but only 0.009 electron plasma periods and negligible proton orbital evolution. It checks execution and short-time numerical behavior; it cannot establish long-term confinement, a sheath equilibrium, beam slowing, or an optimized injection angle.

At 1e10 m^-3, physical Coulomb relaxation is far slower than this horizon. A nearly unchanged velocity distribution is expected even with the collision operators active. Fast visible redistribution must not be manufactured by treating unresolved high-density settings as validated kinetic physics.

## Boundaries and what the diagnostics can establish

This package explicitly selects **standard periodic z boundaries** as an experimental variant while Lars's custom `Periodic_ReflectParVel` implementation is unresolved. Standard periodic wrapping shifts position to the opposite face and leaves the velocity vector unchanged. It is not equivalent to a boundary that rotates velocity, reverses a component, or changes speed. Rotating or reversing a velocity can preserve its magnitude; changing its magnitude also changes kinetic energy. The exact custom map must come from the source that implements it.

The external field's transverse components do not match across the longitudinal seam. Therefore the short kinetic test keeps particles away from that seam; it does not validate trajectories that cross it. Standard periodic wrapping also removes longitudinal end escape as a loss mechanism. Transverse faces and the embedded wall absorb particles, so recorded losses describe those surfaces, not an open-ended magnetic-mirror loss cone.

Use particle counts, weighted particle energies by species, electric-field energy, field extrema, and velocity-space diagnostics together. Continuous injection adds particles and energy; absorbing surfaces remove both. The fixed ion bath constrains orbital motion, and analytic stopping transfers energy to an unrecorded reservoir. Consequently total energy of the particles remaining in the box is not expected to be constant. In the kinetic case an energy balance must include injected and escaped kinetic energy as well as particle and electrostatic-field energy. Prescribed static magnetic fields do no direct work through \(q\mathbf v\cdot(\mathbf v\times\mathbf B)=0\).

The half-timestep run can expose gross timestep sensitivity at a fixed horizon. It does not replace density, grid, particle-count, collision-rate, and longer-time checks. Stochastic runs should be compared through ensemble observables rather than requiring individual trajectories to match.

## Implementation references

The following were checked against local WarpX source commit `bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc`; the executable's actual build revision must be recorded separately by the runner.

- [WarpX input parameters: particles, collisions, electrostatic fields, diagnostics](https://warpx.readthedocs.io/en/latest/usage/parameters.html).
- [Versioned collision wrapper and Coulomb-log handling](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Particles/Collision/BinaryCollision/Coulomb/PairWiseCoulombCollisionFunc.H).
- [Versioned Perez binary-collision implementation](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Particles/Collision/BinaryCollision/Coulomb/ElasticCollisionPerez.H).
- [Versioned background-stopping implementation and kelvin conversion](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Particles/Collision/BackgroundStopping/BackgroundStopping.cpp).
- [Versioned deterministic Cartesian particle positions](https://github.com/BLAST-WarpX/warpx/blob/bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc/Source/Initialization/InjectorPosition.H).

This guide explains the model choices; generated inputs, run logs, and measured diagnostics determine what actually ran.
