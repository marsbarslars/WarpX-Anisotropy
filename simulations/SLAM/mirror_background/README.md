# Continuous mirror injection with plasma backgrounds

This is Ryan's 2026-09-10 extension of Lars's continuous, already-ionized proton
source. It supports the future SLAM work, but the geometry and field here are a
tutorial magnetic mirror. Read `RESULTS.md` for the actual completed runs.

## Start here

1. `workflow.py`: generates a readable WarpX input, launches the existing GPU
   executable, and checks particle/field diagnostics.
2. `PHYSICS.md`: background particles, collisions, collective fields, and math.
3. `ASSETS.md`: located magnetic field, units, checksums, and wall geometry.
4. `BOUNDARY_STATUS.md`: exactly what standard periodic wrapping does, and what
   is still missing from Lars's custom boundary implementation.
5. `runs/<case>_<tag>/inputs.txt`: the complete input that actually ran. Open this
   before reading implementation details. `warpx_used_inputs` is WarpX's parsed
   parameter dump; `metadata.json` records the chosen physical model and assets.

`report.py` refreshes the compact measured `RESULTS.md` and `RESULTS.json` from
the four named startup runs. `test_workflow.py` checks model composition and
units without launching the GPU solver.

Lars's original input is `../../Mirror/inputs_mirror.txt`. The workflow reads its
active proton-source parameters. The old commented Gaussian beam and `Np=1000`
do not control continuous injection.

## Available models

| `--case` | Background | Collisions | Deposited electric field |
| --- | --- | --- | --- |
| `ion_background` | Spatially fixed thermal ions | None | Off |
| `ion_coulomb` | Spatially fixed thermal ions | Beam-ion Coulomb | Off |
| `hybrid` | Fixed-position ions and analytic electron bath | Beam-ion Coulomb and electron drag | Off |
| `kinetic` | Moving thermal ions and physical-mass electrons | All six self/cross species pairs | On |
| `kinetic_half_dt` | Same as kinetic | Same as kinetic | On; half timestep, twice the steps |

The first variant is available for inspection, but another collisionless baseline
is not a prerequisite and was not the task for this session. The two reservoir
collision cases use `n=1e16 m^-3`; the kinetic cases use `n=1e10 m^-3` so the
existing grid resolves the initial screening scale. Their output must not be
compared as if density and physics were matched.

## Run in WSL Ubuntu

The existing Python environment and GPU executable are used directly; no CUDA
rebuild or conda activation is required for the commands below.

```bash
cd /mnt/c/Users/ryanv/Desktop/BOOGERAIDS/WIPPL/WarpX-Anisotropy/simulations/SLAM/mirror_background
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python workflow.py run --case kinetic --boundary standard-periodic
```

Replace the case with `ion_coulomb`, `hybrid`, or `kinetic_half_dt`. By default,
kinetic is 200 steps of 5e-11 s (10 ns); its half-timestep companion is 400 steps
of 2.5e-11 s for the same physical horizon. Reservoir cases are 200 steps of
1e-7 s (20 microseconds). `--steps` overrides the nominal step count; the
half-timestep case doubles that number. `prepare` writes a complete input and
metadata without running. `analyze --run-dir runs/<name>` rechecks saved data.
Run names are timestamped; existing runs are never overwritten.

Defaults may be overridden with `--field /absolute/path/to/field.h5` and
`--warpx /absolute/path/to/warpx.3d`. The field must pass the implemented mesh
coverage/unit checks. For a different asset, also inspect its periodicity,
coordinate system, and vessel alignment; matching extents is not enough.

The executed geometry is an **analytic 0.9-m-radius grounded cylindrical bore**,
matching the central inner-wall radius of Lars's STL approximately. The exact
faceted annular STL failed the self-consistent electrostatic solve at the current
resolution; `--wall stl` remains available to reproduce that issue. The analytic
case removes the exterior fluid and uses a smooth cylinder. It does not silently
claim validation of the original STL discretization.

`--boundary standard-periodic` is explicit because the installed solver does not
implement `Periodic_ReflectParVel`. The available field does not join smoothly at
the z seam. Backgrounds start at r<0.75 m, 0.25<z<4.75 m, and these short tests
check distance from the seam. Longer wrapped trajectories need Lars's source
patch and a compatible field mapping.

## Outputs and interpretation

Every run includes:

- `run.log`, `execution.json`: actual GPU startup, exit code, wall time, binary
  hash, and source checkout commit.
- `inputs.txt`, `warpx_used_inputs`, `metadata.json`: exact parameters/provenance.
- `diags/reduced/`: counts and physical weights, weighted particle energies and
  momenta, field energy, and field maxima at every step.
- `diags/diag/openpmd/`: five full openPMD HDF5 snapshots by default, including
  particle positions/momenta/weights, E, B, charge density, potential, and EB mask.
  These `.h5` files are simulation outputs; the separate field-map `.h5` is input.
- `summary.json`: operational checks and measured observables, including explicit
  limitations. A completed run is not automatically a physics-validation result.

The kinetic configuration includes collective electrostatic fields and Coulomb
collisions. It does not solve for self-consistent magnetic perturbations, simulate
neutral ionization, or demonstrate full collision relaxation during 10 ns.
The inherited beam flux is only 1e5 physical protons/(m^2 s), so its feedback on
the background is tiny. Most short-time E evolution comes from the background's
thermal motion and finite spatial profile, not a strong injected beam.

Large outputs stay in the locally ignored `runs/` directory; compact reports and
source code can be committed separately. Nothing in this workflow pushes Git or
changes Lars's WarpX source.

## Temperature at the mirror center

`TEMPERATURE.md` explains the measurement, interpretation, and what to send Lars.
To collect the central particles at **every timestep**, without writing the full
domain every timestep:

```bash
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python workflow.py run --case kinetic --boundary standard-periodic --temperature-every 1
```

This adds a particle-only diagnostic in `r<0.25 m, 2.25<z<2.75 m`. It changes
output settings only, not injection, collisions, timesteps, or boundaries.
The ordinary full-domain snapshots are retained for independent cross-checks.
`--temperature-every 0` (default) leaves the probe off; a positive integer sets
its cadence. For a longer future run, choose a cadence appropriate to the physics
and disk budget. The probe is spatially filtered but still writes raw particles.

After the run, `temperature.py` automatically produces:

- `temperature_history.png`: ion/electron temperatures and local particle counts.
- `temperature_history.csv`: step, physical timestamp, species, temperature,
  Cartesian component temperatures, mean velocity, bulk energy, and weights.
- `temperature_summary.json`: compact results and comparisons with independently
  selected particles in the ordinary snapshots.
- `diags/center/openpmd/`: raw central particles, all local markers with physical
  weights; no random downsampling. Empty beam selections have blank temperature,
  not zero temperature.

To regenerate these analysis products without rerunning the simulation:

```bash
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python temperature.py runs/kinetic_temperature200_20260914
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python -m unittest -v test_workflow test_temperature
```
