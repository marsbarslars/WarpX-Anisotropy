# Verified hybrid startup results — 2026-09-18

## Latest: evolving electron energy and 1-microsecond diagnostics

Read `TODAY_HANDOFF.md` first for the current completed work, actual input reading
guide, numerical results, and six slide-ready PNG/SVG figure pairs. The latest
run is `runs/energy_transport_1us_20260918`, with a same-horizon half-timestep
reference in `runs/energy_transport_halfdt_1us_20260918`.

Both GPU runs completed: 200×5 ns (126.76 s wall time) and 400×2.5 ns (134.42 s).
Electron entropy/internal energy now evolves, gamma=5/3, with transport and
compression/expansion. No electron-ion or Joule heating source is enabled.
Maximum background-ion timestep difference is 0.033067% of stated initial
scales; electron-fluid difference is 0.000464594%. All 63 hybrid unit tests pass.
This is numerical startup agreement, not steady-state or beam-heating validation.

`VERIFICATION_RESULTS.md` records the isolated collision-operator tests:
thermal relaxation agrees closely, but opposing ion streams lose 1.8651% of
ion-plus-electron energy in the tested relaxation path, without corresponding
electron heating; the zero-rate control changes by only 1.33e-7 percent.
Consequently that installed operator stays off in the mirror.

The sections below retain the earlier fixed-temperature results for provenance.

Scope: true fluid-electron/kinetic-ion solver and diagnostic **operational
validation**, in the explicitly labelled analytic periodic mirror fixture.
Not SLAM vessel results, a heating measurement, or a confinement study.

## Completed GPU runs

| Run folder under `runs/` | Closure | Ion collisions | Steps / dt | Physical time | Measured executable wall time |
| --- | --- | --- | --- | --- | --- |
| `isothermal_20ns_20260918` | Fixed Te=10 eV | Off | 20 / 1 ns | 20 ns | 19.36 s |
| `isothermal_halfdt_20ns_20260918` | Fixed Te=10 eV | Off | 40 / 0.5 ns | 20 ns | 21.55 s |
| `isothermal_ioncoll_50ns_20260918` | Fixed Te=10 eV | On | 50 / 1 ns | 50 ns | 26.81 s |
| `polytropic_4ns_20260918` | gamma=5/3 | Off | 4 / 1 ns | 4 ns | 4.64 s |

All exited with status 0 and reached their requested final output time. Runtime
logs identify **Hybrid-PIC (Ohm's law)** on the collocated grid. There are only
`protons` and `background_ions` particle species. The collisional run log records
50 `doCollisions` calls and both configured binary operators; it is not merely
a label. None of these cases contains electron-ion collisional energy exchange.

The installed executable is
`/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d` (CUDA, double precision,
openPMD, EB); every run manifest records its SHA256 and the source-input SHA256.
No WarpX source or build was changed.

## What the latest 50-ns ion-collision run says

| Final measured quantity | Center | Near injection |
| --- | --- | --- |
| Background markers | 11,232 | 1,728 |
| Background n (whole-cell counts) | 1.0e16 m^-3 | 1.0e16 m^-3 |
| Ion T_parallel | 0.958270 eV | 0.930312 eV |
| Ion T_perp | 0.958558 eV | 0.951671 eV |
| Ion T_scalar | 0.958462 eV | 0.944551 eV |
| Ion P_parallel | 1.535317e-3 Pa | 1.490525e-3 Pa |
| Ion P_perp | 1.535778e-3 Pa | 1.524744e-3 Pa |
| Electron-fluid Te (imposed) | 10 eV | 10 eV |
| Electron-fluid Pe | 0.01602219 Pa | 0.01601921 Pa |
| Beam markers | 0 | 1,034 |

Initial measured scalar ion temperatures were 0.958465 eV and 0.944566 eV,
respectively. The difference from the input 1 eV mostly reflects finite-marker
sampling and per-cell sample-mean subtraction, not time evolution. With 27
markers/cell, the empirical variance has a 26/27 expectation factor for independent
Maxwellian samples. Read the README before interpreting absolute small offsets.

The total background population remains 974,400 markers. The 1,034 beam markers
represent only **0.000201953 physical protons** at this artificially weak source
flux. Macro-particle counts are not physical particle counts or density. No beam
has reached the center. A temperature rise cannot be claimed from these data.

The finite plasma column is not initialized in equilibrium. Large electric
fields near the radial edge are allowed startup responses to the sharp pressure
gradient (the 20-ns test has global max |E| about 228 V/m), not evidence of beam
heating. Isothermal electrons are an imposed thermal reservoir, so total closed
particle+electron energy conservation is not a valid claim for this closure.

## Checks that passed, and what they do not prove

- 20 unit tests: weighted moments, local drift removal, rotated local B,
  anisotropic moments, empty/low-count cells, normalization/tails, mask volumes,
  invalid data, source/model selection, periodic/divergence-free analytic field.
- Completed-run analysis checks finite fields/particles, both ion species only,
  whole-cell coordinate alignment, no EB-covered selected cells, final time,
  constant background count, and captured histogram weight.
- After initialization, `Pe = rho_ion * Te[eV]` holds to roughly 1e-15 relative
  to n0*e*Te in the isothermal selected regions. Native `Te` output was correctly
  converted from Kelvin to eV.
- Halving dt at the same 20-ns endpoint changes selected background temperature
  components by at most **4.46e-8 relative**, and scalar T by at most 1.15e-8.
  This is a short-startup bulk-moment check, not convergence of fast waves,
  heating, collision rates, losses, or the beam PDF. NFlux injects twice as many
  lower-weight markers at half dt, so its random realization differs.
- The polytropic case runs and yields a density-dependent fluid temperature;
  final central Te=10.0000141 eV and near-source Te=9.9999143 eV. This is an
  algebraic pressure response, not electron heating. Unlike gamma=1, averaging
  nodal nonlinear closure quantities to cell centers does not preserve an exact
  pointwise closure equality; no unsupported equality check is claimed.

There has NOT been a density-floor sweep, grid/particle-number convergence,
long-time equilibrium check, electron-ion energy-transfer test, beam-off control,
or real SLAM boundary validation. The normal periodic wrap is not an end-loss
diagnostic. The scalar electron closure does not measure electron anisotropy.

`validate_runs.py` re-creates figures and writes the machine-readable combined
summary to ignored `runs/validation_summary.json`.

## Graphs and numerical records

Use the `isothermal_ioncoll_50ns_20260918` folder for the latest ion-collision
test. All paths below are relative to that run:

- `ion_fluid_history.png`: measured background T/P and imposed fluid Te.
- `velocity_distributions.png`: physical-weighted ion PDFs, separately by region
  and population. The empty central beam panel is intentional.
- `density_and_counts.png`: density, computational markers, represented physical
  particles; symlog axes include empty populations.
- `sampling_cells.png` and `sampling_cells.json`: cell map and complete i/j/k list.
- `ion_history.csv`, `electron_fluid_history.csv`, `velocity_distributions.npz`:
  numerical data. There is no smoothing or invented intermediate sampling.
- `inputs.txt`, `manifest.json`, `warpx_used_inputs`, `warpx.log`, `execution.json`:
  actual settings, identity, log and completion evidence.

## Lars's update found — not merged or executed

Read-only remote branch discovery found `main` and `RyanZhu`. `origin/main` was
fetched (Git metadata only), and the new commit is:

[50d81ac — Included SLAM mirror and custom periodic BCs](https://github.com/marsbarslars/WarpX-Anisotropy/commit/50d81ac9d0bea76b32a4c30ac7d989bab3339ac7),
authored by Lars on September 18, 2026.

It includes the updated `simulations/Mirror/SLAM_VV.stl`, also changes
`simulations/SLAM/SLAM_VV.stl`, and changes the Mirror input substantially:

- Box: x=-0.77..0.77 m, y=0.25..0.75 m, z=-0.25..0.25 m; mesh 96x32x32.
- New source plane y=0.6 m, normal -y; source center x=0.5,z=0, half-widths
  0.05 m; ux,uy=-0.0004..-0.00035; uz=+/-9.5e-5; 4 markers/cell/step.
- All field boundaries PEC; particle z boundary `Periodic_ReflectParVel`.
- Applied field file is now `SLAM_field.h5`, still particle-only in that deck.
- It ships `warpx_modified.3d`, but the commit does not include the corresponding
  WarpX boundary source patch/build recipe. The downloaded executable was NOT run.

`SLAM_field.h5` is not tracked in that commit and was not found by filename in
the local WarpX-Anisotropy, hackathon, or lab-data folders checked. Ask Lars for
its location plus the boundary source branch/patch and exact wrapping rule.
The old center z=2.5 m is **outside the new box**, so reusing its cell selections
would be wrong. A vessel STL alone is not enough to define the new diagnostics
or a consistent fluid-electron field setup.

No merge/pull into the working files, asset overwrite, source rebuild, commit,
or push was performed. Existing work stays on `RyanZhu`. A future merge will
change the live Mirror source, so this fixture generator deliberately rejects
the new source geometry until that migration is reviewed.
