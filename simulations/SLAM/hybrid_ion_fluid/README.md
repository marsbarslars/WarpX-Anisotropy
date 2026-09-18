# Kinetic ions + fluid electrons: separate hybrid-PIC validation

For today's completed-work summary, plots and input reading order, start with
`TODAY_HANDOFF.md`. The older fixed-temperature startup results are retained
in `RESULTS.md` for provenance, not replaced.

This is the new **true WarpX hybrid-PIC** workflow. It does not replace or
rename the five legacy presets in `../mirror_background/`. In particular, that
folder's `hybrid` preset means prescribed fields plus approximate electron
stopping; it is NOT the fluid model here.

## What is actually simulated

- `background_ions`: moving, depositing proton macro-particles, initially 1 eV.
- `protons`: Lars's live continuous, already-ionized proton source, copied from
  `../../Mirror/inputs_mirror.txt`. Beam and background remain separate species.
- Electrons: inertialess, quasi-neutral fluid; no electron particle species.
  Its density follows the total ion charge density, and its pressure/current
  enter Ohm's law. The magnetic perturbation evolves with Faraday's law.
- Initial background density is **1e16 m^-3**, matching the old reservoir
  density, NOT the old kinetic-electron numerical test's 1e10 m^-3.
- The default validation is collisionless (collective fields are still active).
  `--ion-collisions` adds beam/background and background/background Coulomb
  collisions, with ln Lambda = 10. Electron-ion collisional drag/energy exchange
  is NOT enabled by either choice. Extreme macro-weight ratios and the disabled
  collision energy correction must be revisited for quantitative heating.

Small electron mass creates short kinetic electron timescales. A fluid model
removes those resolved electron orbits; it does NOT freeze electrons in place
and does NOT remove electron pressure or electric fields. It also does not
resolve electron velocity distributions, electron Debye sheaths, light waves,
or arbitrary high-frequency electron physics. Hall/whistler field substeps
still restrict the numerics. This installed hybrid solver has one mesh level,
not AMR support.

## Electron closure is a visible decision

| Choice | What changes | Meaning/limitation |
| --- | --- | --- |
| `--closure isothermal` (default) | Te = 10 eV; Pe = ne * e * Te | Te is imposed; no electron-heating result or closed electron energy budget |
| `--closure polytropic` | gamma=5/3; Te = 10*(ne/n0)^(2/3) | Density-dependent algebraic response, NOT evolved heating/transport |
| `--closure energy_transport` | Evolve electron entropy/internal energy, gamma=5/3 | Transport and compression/expansion; no collisional/Joule heat source, conduction, or physical electron wall-loss model |

The evolving-energy option is now implemented. Read `ELECTRON_MODEL_AUDIT.md`
and `VERIFICATION_RESULTS.md` before enabling electron-ion exchange: the
installed relaxation operator does not return independent ion drift energy
to the electron fluid. It stays off in the mirror case. Enabling a fluid alone
does not select a physically justified drag rate. A reconstructed Maxwellian
would be an assumption, not a measured
electron velocity distribution, so we do not plot one.

## This mirror is an explicitly labelled test fixture

The old FEMM B map is applied only to particles in the legacy workflow, and its
off-axis field does not match across the periodic z seam. Simply enabling a
fluid solver with that particle-only field would be inconsistent. Instead,
this separate validation uses a smooth periodic analytic vector potential:

```
A = 0.5 * Bavg * [1 + a*cos(k*z)] * (-y, x, 0)
B = curl(A)
k = 2*pi/(5 m)
Bavg = Bmin*(Rm+1)/2; a = (Rm-1)/(Rm+1)
Bmin = 0.001487675 T; Rm = 3.10951921
```

Only the on-axis field scale and ratio are matched; this is NOT a reconstruction
of FEMM or a SLAM result. The vector-potential field splitting supplies the
applied field consistently to ions and fluid. No particle-only B is added.
The validation bore is analytic r=0.9 m, not Lars's new vessel. Fields use PEC
side boundaries; particles are absorbed at the wall and ordinary-periodic in z.
Recirculation is NOT physical escape through mirror ends or Lars's custom wrap.

Background initially fills r<0.75 m over the full periodic z extent. It is
**not a force-balanced equilibrium**; edge expansion and initialization changes
cannot be attributed to beam heating. The density floor is 1e15 m^-3 and the
Holmstrom low-density treatment is on. These are numerical regularizations,
not added physical background particles, and need sensitivity checks later.

## Exactly which cells are measured

Mesh: 32 x 32 x 80 on x,y in [-1,1] m and z in [0,5] m.
All cell widths are 0.0625 m. Indices are **zero-based**, `(i,j,k)=(x,y,z)`,
and the following upper index bounds are exclusive.

| Region | Selected cells | Actual volume | Why |
| --- | --- | --- | --- |
| Center | i,j in 12:20, retaining cell centers with x^2+y^2<0.25^2; k in 36:44 | 416 cells; 0.1015625 m^3 | Interior, low-field bulk; z=2.25 to 2.75 m, away from source |
| Near injection | i in 12:16, j in 14:18, k in 14:18 | 64 cells; 0.015625 m^3 | x=[-0.25,0), y=[-0.125,0.125), z=[0.875,1.125), includes source and immediate downstream space |

The source plane is x=-0.2 m with y +/-0.1 m and z=0.9..1.1 m.
We include whole cells, NOT arbitrary fractions of cells. The central mask is
a staircase approximation to a cylinder; we use its actual summed cell volume,
not pi*r^2*length. It is intentionally not identical to the old particle-wise
cylindrical filter. `sampling_cells.json` enumerates every selected cell.

Both selections are well inside the bore; analysis refuses cells touching the
embedded boundary. On a new vessel/mesh the masks must be redesigned rather
than silently reusing these indices. Cell width, local B, source location,
wall/field geometry, sufficient markers, and distance from the density-floor
region determine whether a proposed selection is sensible.

## Diagnostic definitions

For each species and each selected cell, let w be physical macro-particle
weight, V the cell volume, and U = sum(w*v)/sum(w). Set c=v-U. Project c onto
the **local total B at that cell center**, not a presumed global z direction:

```
n = sum(w)/V
P_parallel = m_p/V * sum(w*c_parallel^2)
P_perp = m_p/(2*V) * sum(w*(|c|^2-c_parallel^2))
T_parallel[eV] = P_parallel/(n*e)
T_perp[eV] = P_perp/(n*e)
T_scalar = (T_parallel + 2*T_perp)/3
P_scalar = (P_parallel + 2*P_perp)/3
```

Region pressures are volume averages of the cell pressures, and region
temperatures are therefore physical-particle-weighted cell temperatures.
Subtracting one mean velocity for the entire region instead would include
resolved shear in the thermal measurement. Sub-cell flow variations remain
unresolved in this implementation.

These are empirical second moments, NOT Maxwellian fits and not unbiased
finite-marker estimators. With 27 equal-weight markers/cell, removal of the
sample cell mean gives an expected 26/27 factor for ideal independent Maxwellian
samples: an observed value near 0.96 eV can occur even when loading 1 eV.
That does NOT mean physical cooling. We do not secretly apply a Bessel correction
to the pressure. Record marker count and effective count `(sum w)^2/sum(w^2)`;
cells below effective count 10 are flagged. More markers/coarser measurement
bins and convergence checks are needed for small heating signals.

Velocity histograms show **absolute lab-frame** v_parallel and |v_perp| using
the same local B. Beam and background are normalized separately by physical
weight. This is a radial v_perp probability density, not Cartesian f(vx,vy,vz);
no division by 2*pi*v_perp is applied. Bin limits are +/-100 km/s parallel and
0..100 km/s perpendicular. The captured weight fraction is recorded so tails
outside the plot cannot disappear unnoticed. Empty populations have undefined
temperature (blank), not a zero-temperature claim.

Fluid output `Te` is in **Kelvin in WarpX** and converted to eV in analysis;
`Pe` is Pa. Grid electron density follows deposited/smoothed ion charge, so it
need not equal a nearest-cell ion-count estimate exactly. Iteration 0 is dumped
before hybrid field/pressure initialization. At that iteration only, initial
particle projections use the analytic applied B; initial fluid pressure/density
are flagged uninitialized, not interpreted as vacuum. Later projections use
the saved total simulated B without adding the applied field a second time.

For an evolving fluid, the CSV reports volume-mean Te, its cell minimum/maximum,
and `Te_effective_eV = <Pe>/(e*<ne>)`. These means need not be identical; the
nodal fields are averaged to cells independently. Regional electron internal
energy is `sum(Pe*cell_volume)/(5/3-1)`, an inventory rather than a complete
open-mirror conservation budget. Below-density-floor cell fractions are saved.
The saved fluid velocity supplies a cell-centered advection-CFL proxy, not a
strict bound on internal nodal or entropy-marker motion.

## Files and usage

- `setup_case.py`: commented input generation and opt-in GPU execution.
- `diagnostics.py`: moments, exact cell masks, weighted distributions, plots.
- `test_diagnostics.py`: numerical, indexing, closure-configuration tests.
- `test_setup_case.py`: explicit energy/cadence/adaptive-step configuration tests.
- `presentation_plots.py`: readable PNG/SVG figures and matched-timestep differences.
- `verify_energy.py`: isolated closed periodic-box thermal/drift-energy tests;
  artificial relaxation rate, not a prediction for the mirror.
- `ELECTRON_MODEL_AUDIT.md`: installed-version equations and limitations.
- `VERIFICATION_RESULTS.md`: measured isolated energy-operator checks.
- `RESULTS.md`: verified runs and caveats; not generated promotional claims.
- `runs/<name>/inputs.txt`: actual full, readable WarpX input for that run.
- `runs/<name>/ion_history.csv`: both species, both regions, T/P/n/counts.
- `runs/<name>/electron_fluid_history.csv`: Te, Pe, grid density/closure checks.
- `runs/<name>/velocity_distributions.npz`: histogram arrays, shared bin edges.
- `runs/<name>/sampling_cells.json`, `diagnostic_audit.json`: selections/checks.
- PNGs in each run: history, distributions, selected-cell visualization.

Use the existing WSL `warpx-gpu` environment; no rebuild or activation needed
when using its interpreter and binary paths. From this folder in WSL:

```bash
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python -m unittest -v test_diagnostics
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python setup_case.py runs/my_isothermal --execute
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python diagnostics.py runs/my_isothermal
```

Add `--ion-collisions`, `--closure polytropic`, or `--closure energy_transport`
explicitly. Existing run folders are refused. Default dt is 1 ns; up to 5 ns
requires `--adaptive-fields`, which adapts magnetic-field substeps only (not
particle or fluid-advection steps). This bound is a workflow guard, not a
universal stability guarantee. Check saved fluid CFL proxy and half-dt results.
Use `--diagnostic-every N` to specify saved-output cadence independently of
the time horizon. Fixed field subcycling is requested at 64; adaptive RKF45 uses
rtol=1e-5, atol=1e-12 and records accepted/rejected substeps in the runtime log.

Commit code, tests, and reviewed small notes; raw runs remain ignored. Do not
commit WarpX source/builds, large HDF5 diagnostics, or videos. Nothing here is
pushed automatically.

## Before moving onto Lars's vessel

1. Confirm electron closure and whether electron-ion drag/heating is required.
2. Obtain the exact remote branch/commit for the fixed vessel and compatible
   boundary implementation, plus consistent applied field data (preferably A).
3. Inspect changes without overwriting the existing dirty working tree.
4. Redesign whole-cell/volume masks against the actual mesh and vessel; use
   local B for parallel/perpendicular moments and handle cut cells explicitly.
5. Define a useful beam flux/energy, run duration, equilibrium/background, and
   a beam-off comparison. Lars's inherited 4000 physical protons/s is only an
   approximately 8e-15 W source here: the current test is not a heating study.

Primary reference: [WarpX kinetic-fluid hybrid model](https://warpx.readthedocs.io/en/latest/theory/models_algorithms/kinetic_fluid_hybrid_model.html).
Runtime input syntax was also checked against the installed WarpX source.
