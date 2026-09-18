# SLAM WarpX study guide for VS Code

This file is a reading guide for Ryan, Lars, Tony, and a future VS Code coding
assistant. It describes what is currently implemented, which files to open, what
symbols to inspect, and which physics questions to answer before changing the
simulation. It is **not** permission to launch a new vessel run, overwrite an
input, replace an STL, or push to Git.

## The one-sentence project context

We are studying a neutral-beam-like proton population moving through a magnetic
mirror while a kinetic-ion / inertialess-electron hybrid WarpX solver evolves
the background plasma. The current case is a controlled analytic mirror fixture;
it is not yet the SLAM vessel case.

## What is complete

- The existing analytic periodic mirror is retained: a smooth vector potential
  creates a periodic magnetic field with on-axis mirror ratio `Rm = 3.10951921`.
- `background_ions` are mobile kinetic protons, initially nominally 1 eV and
  `n0 = 1e16 m^-3`, filling `r < 0.75 m`.
- `protons` are a separate continuous injected beam copied from the old mirror
  source: about 12 eV mean energy, flux `1e5 m^-2 s^-1`, about 4000 represented
  physical protons per second.
- Electrons are a **fluid**, not particles. The current `energy_transport`
  closure evolves electron internal energy with `gamma = 5/3`; it includes
  transport and compression/expansion, but no heat conduction, kinetic-electron
  velocity distribution, or physical electron sheath/wall-loss model.
- Electron-ion collisional drag/heating is intentionally **off**. The installed
  relaxation operator failed a directed-energy accounting test: counterstreaming
  ions lost drift energy without a corresponding electron-energy increase.
- The final runs reached 1 microsecond at 5 ns and 2.5 ns timesteps. Their
  matched-sample differences are numerical checks, not experimental accuracy,
  steady-state evidence, or proof of beam heating.
- The local tests passed: 63 hybrid/model/plot tests and 9 earlier temperature
  tests.

Read `TODAY_HANDOFF.md` for the measured values and plain-language explanation.
Read `ELECTRON_MODEL_AUDIT.md` for the installed-source equations and the
collision-operator warning. Read `VERIFICATION_RESULTS.md` for the isolated
thermal and counterstream tests.

## Open these files in this order

Use the Windows workspace root:

`C:\Users\ryanv\Desktop\BOOGERAIDS\WIPPL\WarpX-Anisotropy`

1. `simulations/SLAM/hybrid_ion_fluid/VSCODE_STUDY_GUIDE.md` — this map.
2. `simulations/SLAM/hybrid_ion_fluid/TODAY_HANDOFF.md` — what was run today,
   initial conditions, diagnostic regions, results, and limitations.
3. `simulations/SLAM/hybrid_ion_fluid/runs/energy_transport_1us_20260918/inputs.txt`
   — the exact WarpX input that ran on the GPU.
4. `simulations/SLAM/hybrid_ion_fluid/setup_case.py` — how that input is built.
5. `simulations/SLAM/hybrid_ion_fluid/diagnostics.py` — how every plotted
   temperature, pressure, density, PDF, and count is calculated.
6. `simulations/SLAM/hybrid_ion_fluid/presentation_plots.py` — how the
   slide-ready PNG/SVG figures and timestep misfits are made.
7. `simulations/SLAM/hybrid_ion_fluid/runs/energy_transport_1us_20260918/`
   — CSV/JSON data, runtime log, execution status, and `presentation/` figures.
8. `simulations/SLAM/hybrid_ion_fluid/RESULTS.md` — earlier fixed-temperature
   results plus the current short summary.
9. `simulations/SLAM/hybrid_ion_fluid/ELECTRON_MODEL_AUDIT.md` and
   `VERIFICATION_RESULTS.md` — model-validity reading, not promotional claims.

Do not confuse the new workflow with the earlier hackathon workflow in
`simulations/SLAM/mirror_background/`. That folder is useful historical code,
but its old `hybrid` preset is not this true fluid-electron solver.

## Project map

| Path | Role | What to learn |
|---|---|---|
| `simulations/SLAM/PROJECT_CONTEXT.md` | Project-level handoff | Why the work moved from hackathon scans to SLAM, and what is still provisional |
| `simulations/SLAM/hybrid_ion_fluid/setup_case.py` | Input generator | Domain, field, species, closure, source, diagnostics, manifest |
| `simulations/SLAM/hybrid_ion_fluid/diagnostics.py` | Analysis | Cell masks, local-B moments, fluid fields, PDFs, inventories, CFL proxy |
| `simulations/SLAM/hybrid_ion_fluid/presentation_plots.py` | Figure generator | Temperature/pressure/count/PDF/misfit figures without interpolation |
| `simulations/SLAM/hybrid_ion_fluid/verify_energy.py` | Isolated verification | Artificial-rate thermal test and counterstream drift-energy test |
| `simulations/SLAM/hybrid_ion_fluid/test_*.py` | Regression tests | What the implementation promises mathematically and structurally |
| `simulations/SLAM/inputs/inputs_3d_magnetic_mirror.txt` | Older tutorial-style input | Historical test-particle mirror; not the current fluid case |
| `simulations/SLAM/inputs_SLAM.txt` | Older SLAM test-particle deck | Not the current hybrid-fluid input; uses a particle-only field path |
| `simulations/SLAM/SLAM_VV.stl` | Local vessel mesh | Current working copy is older than Lars's fetched updated Git object |
| `simulations/SLAM/mirror_background/BOUNDARY_STATUS.md` | Boundary audit | Why ordinary periodic wrapping is not Lars's custom velocity-transform boundary |

## What to inspect in the Python code

### `setup_case.py`

Open these functions and answer what each input line does:

- `build_input(...)`: the source of truth for the generated deck.
- The `energy_equation` / `closure_comment` block: why evolving `Te` is not the
  same as enabling electron-ion heating.
- The applied-vector-potential block: how `A` gives a consistent field to both
  the ion particles and the fluid solver.
- The background-ion initialization: density, Gaussian thermal velocity, and
  marker count.
- The injected `protons` source copied from `simulations/Mirror/inputs_mirror.txt`.
- The diagnostic block and `manifest` dictionary: what metadata is recorded so
  a later comparison cannot silently mix different physics.
- The CLI arguments `--closure`, `--diagnostic-every`, and `--adaptive-fields`.

Questions to write down while reading:

1. Which settings make electrons fluid rather than kinetic particles?
2. Which settings make the background protons kinetic and depositing?
3. Which settings are physical assumptions, and which are numerical safeguards?
4. Which fields are evolving, and which are only initialized or diagnosed?

### `diagnostics.py`

Read these symbols in sequence:

1. `Grid`, `centers`, and `ids`: how a particle is assigned to a Cartesian cell.
2. `regions`: why the center and near-injection ROIs have exact cell lists and
   whole-cell volumes.
3. `moments`: local mean-flow subtraction, projection onto the local total `B`,
   and definitions of `T_parallel`, `T_perp`, `P_parallel`, and `P_perp`.
4. `fluid_moments`: volume-mean `Te`, effective `Te = <Pe>/(e<ne>)`, pressure,
   density-floor flags, and electron internal-energy inventory.
5. `kinetic_energy`: the particle kinetic-energy calculation used for inventory
   checks; it is not by itself a closed energy-conservation proof.
6. `field_xyz`: raw openPMD field names, axes, cell-center checks, and the
   Kelvin-to-eV conversion for WarpX's `Te` output.
7. `analyze`: the iteration loop that reads fields/particles and writes the
   CSV, NPZ, JSON, and PNG outputs.
8. `plot`: which curves are diagnostic summaries and which are not claims of
   equilibrium or heating.

Pay particular attention to these distinctions:

- Ion temperature removes each cell's mean flow; the beam's directed energy is
  not automatically called “temperature.”
- The electron fluid has scalar pressure; it does not provide an electron
  velocity PDF.
- Computational marker count and represented physical-particle count differ by
  macro-particle weight.
- `Te` and `Pe` at iteration zero are excluded from fluid evolution plots because
  the first dump precedes hybrid pressure initialization.
- The saved electron-advection CFL is a cell-centered proxy, not a strict bound.

### `presentation_plots.py`

Inspect `load_run`, `compare_background`, `load_pdf_payload`, and `make_plots`.
Confirm that:

- comparisons use matching saved times only;
- no interpolation or smoothing is introduced;
- the reference run is a timestep reference, not “truth”;
- physical-weight PDFs retain the captured-weight fraction;
- fluid and ion quantities use clear units and separate scales.

## WarpX source code to inspect in WSL

The installed GPU build uses this source tree:

`/home/ryanv/src/warpx`

The executable is:

`/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d`

The audited source commit is `bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc`.

Open these files, preferably with the commit shown in the file header or Git
history:

- `Source/FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.cpp`
  - search `QDSMCInitializeUe`: how the fluid velocity is obtained from current
    and charge density;
  - search `QDSMCInitializeKe`: how the entropy/internal-energy state starts;
  - search `QDSMCUpdateTe`: how the electron temperature is reconstructed;
  - search `QDSMCFillElectronPressureFromTe`: how `Pe` is formed;
  - search `electron_ion_relaxation_rate`, `QDSMCApplyIonHeating`, and the
    relaxation update: why the installed beam-drag path is not yet accepted;
  - search `include_joule_heating`: why it is explicitly zero in the baseline.
- `Source/Fluids/QdsmcParticleContainer.cpp`
  - search `PushX` and the entropy-marker creation/transport; inspect what the
    density floor, periodic direction, and nonperiodic edge actually do.
- `Source/Particles/PhysicalParticleContainer.cpp`
  - search temperature deposition and the local mean-flow subtraction used for
    ion thermal moments.
- `Source/Diagnostics/FullDiagnostics.cpp`
  - search the `Te`, `Pe`, `eb_covered`, and `additional_fields_to_plot` cases;
  - inspect the explicit raw names recorded in the manifest:
    `hybrid_electron_velocity_fp[dir=x]`, `[dir=y]`, `[dir=z]`.

For each function, log four things in plain language: inputs, output field or
particle quantity, units, and which physical equation or approximation it
represents. Do not only copy code; write one sentence about the physics.

## Physics study checklist

### Fluid electrons

- Electrons are not “one giant particle.” They are represented by density,
  velocity/current, pressure, and an entropy/internal-energy variable on the
  mesh.
- The small electron mass motivates removing electron orbits, but it does not
  prove that kinetic tails, sheaths, or heat conduction are irrelevant.
- In this baseline, `Te` can change through fluid transport and compression. That
  is different from saying the beam heated electrons through collisions.

### Ion diagnostics

- `T_parallel` and `T_perp` are second moments around each cell's mean ion flow.
- The parallel direction is the **local total magnetic field**, not always the
  global z-axis.
- The near-source beam can be sparse and nonthermal; low-occupancy flags matter.
- A flat background temperature curve means “little measured change over this
  short run,” not “the plasma is proven steady.”

### Boundaries and energy

- Ordinary periodic z wrapping preserves the full velocity vector. It is not
  Lars's custom boundary that may transform velocity between ends.
- Absorbing radial/embedded boundaries and continuous injection make the mirror
  an open source/loss system. A raw particle-plus-fluid energy inventory is not
  a complete conservation budget.
- Before claiming confinement or heating, specify source power, wall loss,
  field-energy accounting, boundary transformation, and a beam-off control.

## Future STL/SLAM-vessel checkpoint — do not run yet

The local mesh is:

`simulations/SLAM/SLAM_VV.stl`

It is the older working copy. Lars's newer matching Mirror/SLAM STL exists in
fetched Git history from commit `50d81ac`, but it has not been silently copied
over. The SLAM field asset `SLAM_field.h5` and the source patch for the custom
periodic boundary still need to be identified and version-matched.

Before any vessel run, ask Lars for:

1. the exact branch/commit containing the boundary implementation;
2. the matching WarpX/AMReX build and executable;
3. the exact location and schema of `SLAM_field.h5`;
4. whether the field file is `B`, vector potential `A`, current, or particle-only
   data, and its coordinate units/order;
5. the intended velocity transformation at the periodic seam;
6. a new mesh-aware ROI and embedded-boundary treatment.

Do not reuse the current `z=2.25..2.75 m` center ROI on the SLAM mesh. Do not
assume the old `inputs_SLAM.txt` is compatible with the current fluid model.

## What to record in a study log

For each code section, make a short note with this format:

```text
File/function:
Input parameters and units:
Output field/quantity and units:
Equation or physical interpretation:
Numerical approximation:
What this does NOT prove:
Question for Lars/Tony:
```

Minimum quantities to record from every run:

```text
closure, gamma, Te0, Ti0, n0, density floor
grid, ROI bounds, Bmin, mirror ratio, boundary model
dt, total steps, physical horizon, field substep/adaptive settings
beam flux, source area, beam energy, background initialization
Te, Pe, ne, Te_effective, ion T_parallel/T_perp, ion P_parallel/P_perp
marker count, physical count, captured PDF weight, low-occupancy flags
max E, max B, fluid-CFL proxy, exit code, wall time
```

## Copy-paste prompt for a VS Code study session

```text
Read simulations/SLAM/hybrid_ion_fluid/VSCODE_STUDY_GUIDE.md first, then read
TODAY_HANDOFF.md and the exact input at
runs/energy_transport_1us_20260918/inputs.txt.

This is a study-only session. Do not edit files, launch WarpX, replace the STL,
change the field HDF5, commit, or push unless I explicitly ask. Explain the
existing code and physics in plain language and cite the file/function being
discussed.

Walk me through setup_case.py -> diagnostics.py -> presentation_plots.py in that
order. For each important function, tell me: inputs, outputs, units, governing
equation, numerical approximation, and what conclusion is unsafe to draw.

Then inspect the installed WarpX symbols named in the study guide, especially
QDSMCInitializeUe, QDSMCInitializeKe, QDSMCUpdateTe,
QDSMCFillElectronPressureFromTe, the electron-ion relaxation update,
QdsmcParticleContainer PushX, and FullDiagnostics Te/Pe output handling.

Answer these questions explicitly:
1. Why are protons kinetic while electrons are fluid?
2. Why can Te evolve even when electron-ion heating is disabled?
3. How are ion parallel/perpendicular temperature and pressure calculated?
4. Why are marker counts not physical particle counts?
5. Why is the current result not yet a steady-state or heating claim?
6. What evidence is needed before moving to the SLAM STL vessel?

End with a compact study log using the template in the Markdown guide. Keep the
current analytic mirror separate from the future SLAM-vessel work.
```
<!-- Verified against the 2026-09-18 energy-transport handoff. -->

