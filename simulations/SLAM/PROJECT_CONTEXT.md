# SLAM WarpX Optimization - Project Context

## Read this first

This directory is transitioning from a collection of hackathon mirror examples into a semester-long SLAM optimization project using WarpX.

The imported hackathon files are useful as tested patterns for input generation, parameter scans, collision models, diagnostics, analysis, and visualization. They are **not SLAM simulation results** and should not be presented as validation of the SLAM configuration.

This file is the handoff document for future Codex sessions and collaborators. Update it whenever a major assumption, input asset, validation result, or project objective changes.

## Current direction - 2026-09-10 (supersedes the earlier staged plan below)

Ryan asked to build directly on Lars's continuous proton source, add ion
backgrounds/Coulomb collisions/analytic electron drag, and reach a kinetic
electron-plus-ion electrostatic case now, validating while implementing.
Do not repeat Lars's original beam run or a new collisionless tracer as an entry
requirement. The immediate work stays in the mirror as preparation for SLAM.

The implementation is under `simulations/SLAM/mirror_background/`:

- `README.md`: file map and runnable WSL commands.
- `workflow.py`: source reuse, input generation, GPU execution, diagnostics.
- `RESULTS.md`: actual completed runs, numerical checks, remaining limits.
- `PHYSICS.md`: background definitions, collision/field math, units and scales.
- `ASSETS.md`: found legacy field map, hashes, STL geometry, coordinates.
- `BOUNDARY_STATUS.md`: installed source and missing custom boundary evidence.
- `runs/`: ignored local inputs, logs, full openPMD diagnostics and summaries.

Key findings: the existing CUDA binary works; Lars's checked-in Mirror input
enables electrostatic deposition; the field map was found among legacy assets
outside this repository; `Periodic_ReflectParVel` is not in installed upstream
WarpX or the inspected public branches. Its implementation must come from Lars.
The new cases explicitly use standard periodic z, and the field is not smoothly
periodic off-axis. Backgrounds start in an interior column and runs are short.

The STL annular shell caused a Poisson residual plateau in the new kinetic case.
A grounded analytic cylinder with the same 0.9-m central bore radius converges
at the requested strict tolerance. It is an approximation of the faceted inner
wall, with the exterior made solid; it is not validation of Lars's STL solve.

Reservoir collision cases use `n=1e16 m^-3`, `Ti=1 eV`, and analytic `Te=10 eV`
where applicable, with deposition disabled. The kinetic case uses mobile physical
electrons and proton ions, `ne=ni=1e10 m^-3`, `Te=10 eV`, `Ti=1 eV`, and charge
deposition enabled. Its 5e-11-s timestep resolves electron gyromotion. The
10-ns startup and half-timestep companion exercise the implementation, not
long-time confinement or collisional equilibration. No analytic electron drag is
added on top of kinetic beam-electron collisions.

The previous September 4 audit and fidelity ladder below are historical context,
not an instruction to restart the rejected baseline sequence. Use the new
package's measured report before choosing the next run.

## Current repository state

- Repository: `marsbarslars/WarpX-Anisotropy`
- Working branch at the time of this audit: `RyanZhu`
- Branch status at the time of this audit: clean and tracking `origin/RyanZhu`
- Ryan's imported commit: `924d095` (`Add SLAM mirror campaign analysis tools and results`)
- Lars's `main` branch remains the shared baseline.
- Large field maps and raw WarpX outputs are intentionally not stored in Git; they must be shared separately.
- The repository currently has no Python dependency manifest and no root `.gitignore`.
- The root `README.md` currently contains only the repository title.

## Project mission

Use WarpX to study and eventually optimize energetic-particle confinement in a SLAM magnetic configuration.

The immediate scientific workflow should be:

1. Establish a reproducible, collisionless test-particle baseline in the actual SLAM field and vessel geometry.
2. Validate field coordinates, particle orbits, loss classification, numerical convergence, and diagnostics.
3. Scan beam/source parameters such as injection angle, energy, birth position, and spatial width.
4. Add collisional physics in controlled increments and compare against the collisionless baseline.
5. Define an optimization objective and constraints, then run systematic or algorithmic searches.

The optimization target is not yet final. Candidate objectives include maximizing retained fast-ion fraction or retained fast-ion energy, minimizing vessel/end losses, controlling velocity-space anisotropy, or balancing confinement against wall loading.

## Current directory map

```text
simulations/SLAM/
|-- PROJECT_CONTEXT.md              # This handoff and planning document
|-- README.md                       # Brief description of the imported bundle
|-- inputs_SLAM.txt                 # Lars's initial SLAM test-particle deck
|-- SLAM_VV.stl                     # SLAM vacuum-vessel surface
|-- inputs/
|   `-- inputs_3d_magnetic_mirror.txt
|-- scripts/                        # Ryan's imported hackathon campaign scripts
|-- summaries/                      # Compact legacy mirror results and parameters
`-- figures/                        # Selected legacy mirror figures/previews

scripts/                            # Shared repository-level tools from Lars's main
|-- warpx.py                        # Build/run driver; see current limitation below
|-- nbi_scan.py                     # NBI injection-angle scan
|-- racetrack_scan.py               # SLAM racetrack energy/geometric-loss analysis
|-- vessel_scan.py                  # Vessel-radius and wall-loss scan
|-- plot_mirror.py                  # Mirror field, velocity space, and loss-cone checks
|-- animate_particles.py            # PyVista/OpenPMD trajectory visualization
`-- plasma/                         # Shared field, particle, sweep, and chart helpers
```

## Provenance: what was already present and what Ryan added

### Present on Lars's shared baseline

- `simulations/SLAM/inputs_SLAM.txt`
- `simulations/SLAM/SLAM_VV.stl`
- Repository-level WarpX utilities and shared analysis modules under `scripts/`
- Existing mirror simulation and visualization material under `simulations/Mirror/`

The repository-level tools contain the most SLAM-specific work currently available. In particular:

- `scripts/racetrack_scan.py` traces a closed SLAM magnetic axis and studies the geometric gyroradius limit in a racetrack leg.
- `scripts/nbi_scan.py` scans injection angle and compares simulated confinement with a guiding-center loss-cone prediction.
- `scripts/vessel_scan.py` separates axial losses from vessel-wall losses.
- `scripts/plasma/field.py` loads and interpolates field maps and calculates mirror-ratio/loss-cone quantities.
- `scripts/plasma/particles.py` loads OpenPMD particles and derives energy, pitch, and parallel/perpendicular velocity.
- `scripts/animate_particles.py` creates interactive or rendered 3D particle/field visualizations.

### Added on Ryan's branch

Ryan imported the reusable parts of the 2026 hackathon magnetic-mirror campaign into:

- `simulations/SLAM/scripts/`
- `simulations/SLAM/summaries/`
- `simulations/SLAM/figures/`
- `simulations/SLAM/inputs/`

The imported campaign demonstrates:

- parameterized WarpX input generation;
- collisionless, ion-only, electron-stopping-only, and hybrid collision configurations;
- angle and density scans;
- survival, energy, pitch-space, and loss analysis;
- 3D and velocity-space visualization.

These are reusable methods. Their saved numbers and figures describe the tutorial magnetic mirror, not SLAM.

## Existing inputs

### `inputs_SLAM.txt`

This is the closest current starting point for a SLAM run:

- 3D domain: `2 m x 2 m x 5 m`
- Grid: `40 x 40 x 40`
- Timestep: `4.4e-7 s`
- Horizon: 500 steps
- Species: 1,000 protons
- Source: Gaussian spatial distribution centered at `(0, 0, 2.5 m)`
- Momentum: uniform component-wise distribution
- Model: collisionless test particles (`do_not_deposit = 1`)
- Boundaries: absorbing box boundaries
- Field: prescribed from `example-femm-3d.h5`
- Diagnostics: OpenPMD every step

Important limitations:

- `example-femm-3d.h5` is not present in the repository.
- `SLAM_VV.stl` is not referenced as an embedded boundary in this input.
- Losses therefore occur at the rectangular domain boundary, not at the physical SLAM vessel.
- The input does not contain Coulomb collisions, kinetic electrons, or self-consistent collective fields.

### `inputs/inputs_3d_magnetic_mirror.txt`

This is a legacy tutorial-mirror input, not a SLAM input. It follows one proton in a time-scaled prescribed mirror field and is useful only as a numerical/visualization reference.

It also references an absent `example-femm-3d.h5` field map.

## Legacy hackathon benchmark

The imported campaign used a tutorial-scale proton beam and a prescribed magnetic mirror:

- Initial beam energy: about `5.63 eV`
- Typical beam count: 1,000 for validation and 5,000 for campaign cases
- Angles represented in the full saved comparison: 15, 25, 30, 35, 40, 50, 60, and 75 degrees
- Background densities: `1e15`, `1e16`, and `1e17 m^-3`
- Background ion temperature: `1 eV`
- Background electron temperature: `10 eV`
- Coulomb logarithm: 10
- Physics variants: collisionless, pairwise beam-ion Coulomb collisions, analytic electron stopping, and hybrid
- Charge/current deposition disabled, so the model had no self-consistent collective plasma response

Use these results to test analysis and plotting code. Do not use their confinement values as predictions for SLAM.

## Known portability and integration gaps

The files under `simulations/SLAM/scripts/` were reorganized without yet adapting their path assumptions. They should be treated as reference code until these issues are fixed:

1. `generate_inputs.py` expects `parameters.json` beside the script, but the file is currently under `summaries/`.
2. Several scripts expect `runs/`, `visualizations/`, and `example-femm-3d.h5` beside the scripts.
3. `analyze_case.py` and `analyze_campaign.py` try to import `Plasma-Hackathon/scripts`, which does not exist at that calculated path in this repository.
4. `animate_campaign_3d.py` tries to import the old `warpx/mirror_angle_scan/visualization` module, which is not present here.
5. Generated launchers hard-code `/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d`.
6. `scripts/warpx.py` expects a project containing `vendor/warpx` and `pyproject.toml`; neither is currently present in this repository.
7. No raw OpenPMD campaign data is present, so the imported analysis/animation scripts cannot reproduce the saved figures from this clone alone.

The next implementation session should first establish one portable directory convention rather than patching each missing path individually.

## Required external assets

The exact field-map filename and transfer location need to be agreed with Lars. Current code refers to two different names:

- `example-femm-3d.h5` in the existing SLAM and tutorial inputs
- `SLAM_vC5_warpX.h5` in `scripts/racetrack_scan.py`

Before running particles, verify:

- which file is the canonical SLAM magnetic field;
- field units and component ordering;
- grid spacing and global offset;
- coordinate handedness and axis definitions;
- alignment between the HDF5 field grid and `SLAM_VV.stl`;
- whether the field covers the full particle domain;
- a checksum for the external file.

Large HDF5 maps, OpenPMD/BP5 diagnostics, logs, and rendered movies should remain outside Git.

## Physics fidelity ladder

Do not add every physical effect at once. Use this sequence so each change has a meaningful control:

1. **Single-particle field check:** prescribed SLAM magnetic field, no collisions, no deposition.
2. **Collisionless ensemble:** many test ions with a reproducible source distribution.
3. **Physical vessel:** activate or otherwise enforce `SLAM_VV.stl` and classify end versus wall losses.
4. **Source scans:** injection angle, beam energy, birth position, and beam width.
5. **Collisional controls:** ion-only, electron-stopping-only, hybrid, and matched collisionless cases.
6. **Collective fields, only if justified:** kinetic electrons/ions and charge-current deposition require quasineutral initialization and much stricter spatial/time resolution.

The Egedal/FBIS Fokker-Planck paper is conceptual guidance for slowing, pitch-angle scattering, loss-cone transport, and velocity-space distributions. Reproducing its full orbit-averaged operator is a separate modeling path and is not required for the first SLAM test-particle baseline.

## Validation gates

A run should not become part of an optimization campaign until it passes these checks:

1. The field and vessel occupy the same coordinate system and visibly overlap.
2. On-axis/field-line `|B|`, mirror ratio, and expected loss-cone angle are documented.
3. A single collisionless particle conserves kinetic energy to the chosen tolerance.
4. Magnetic moment is acceptably conserved where guiding-center adiabaticity should hold.
5. Halving the timestep does not materially change the confinement/loss result.
6. Increasing particle count does not materially change ensemble statistics.
7. Losses are classified by physical boundary: mirror/end, vessel/wall, or numerical box.
8. Every scan records its input parameters, Git commit, WarpX version/build, random seed, and output location.

## Decisions needed before the first new run

Discuss these with Ryan and Lars before implementation:

1. What is the canonical SLAM HDF5 field map, and how is it distributed?
2. Should `SLAM_VV.stl` be an actual WarpX embedded boundary or initially only a post-processing loss surface?
3. Which energetic species is the baseline: proton, deuteron, or another ion?
4. What beam energy, angle convention, source location, and spatial distribution represent the intended experiment?
5. Is the first milestone purely collisionless test-particle confinement?
6. What scalar quantity will eventually be optimized?
7. What physical and engineering constraints bound the search?

## Recommended next work session

The next coding session should not begin with a large parameter scan. It should:

1. Read this file and inspect `git status` without changing branches unexpectedly.
2. Obtain and identify the canonical SLAM field map outside Git.
3. Choose a portable data/output layout and add environment/dependency documentation.
4. Repair or replace the moved hackathon scripts so they use repository-relative paths and the shared `scripts/plasma/` modules.
5. Create one minimal collisionless SLAM input.
6. Visualize the field together with `SLAM_VV.stl` before launching particles.
7. Run one single-particle smoke test, then a small reproducible ensemble.
8. Record validated baseline parameters and results in this document.

## Starter prompt for a new Codex session

```text
Continue the semester-long WarpX work on branch RyanZhu. Read simulations/SLAM/PROJECT_CONTEXT.md and simulations/SLAM/mirror_background/{README,RESULTS,PHYSICS,ASSETS,BOUNDARY_STATUS}. The immediate project is Lars's continuous mirror injection with ion backgrounds and kinetic electrons, supporting future SLAM optimization. Inspect current Git state and existing validated startup runs. Do not repeat the old collisionless baseline. Preserve the distinction between dense prescribed-background collision cases and the low-density self-consistent kinetic demonstration. Lars's custom periodic velocity mapping and exact field provenance remain unresolved; the analytic cylinder is a documented geometry approximation. Use actual input/log/diagnostic evidence, discuss the intended density/source strength/horizon, and extend the present implementation with appropriate resolution checks.
```

## Update log

- 2026-09-04: Initial context audit created on Ryan's branch. Documented repository provenance, current inputs, reusable tools, missing external field data, moved-script path breakage, validation gates, and the staged path toward SLAM optimization.
- 2026-09-10: Ryan authorized direct background/collision/kinetic-electron work with concurrent validation. Added mirror_background package, local field/source audit, and updated the current direction; see RESULTS.md for execution evidence.
