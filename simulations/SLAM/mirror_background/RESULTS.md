# Measured startup results - 2026-09-10

The four extensions completed on the existing CUDA WarpX build. These results
establish an operational implementation with initial numerical checks. They do
not establish long-time beam confinement, correct collision rates, or a SLAM
optimization result. No new reproduction of Lars's original beam run was needed.

| Run directory under `runs/` | Physical horizon | WarpX wall time | Final beam macroparticles | Mean beam energy (eV) | Status |
| --- | --- | --- | --- | --- | --- |
| `ion_coulomb_smoke200` | 20 microseconds | 13.08 s | 4,135 | 11.723469 | integration_smoke_pass |
| `hybrid_smoke200` | 20 microseconds | 13.07 s | 4,123 | 11.726034 | integration_smoke_pass |
| `kinetic_smoke200` | 10 nanoseconds | 33.07 s | 4,114 | 11.789856 | integration_smoke_pass |
| `kinetic_half_dt_smoke400` | 10 nanoseconds | 60.72 s | 8,255 | 11.762016 | integration_smoke_pass |

Wall times exclude input generation and diagnostic analysis. The machine uses
the existing RTX 5070-family GPU in WSL. Each run log explicitly records CUDA
initialization and completion. The four final runs required approximately
120.0 seconds of solver time in total.

## What ran

- The original continuous source parameters are retained: proton injection from
  x=-0.2 m over y=(-0.1,0.1) m and z=(0.9,1.1) m, physical flux 1e5 m^-2 s^-1,
  component-wise uniform momenta giving approximately 9-16 eV birth energies.
- All final runs use the analytic radius 0.9-m grounded bore and explicitly
  standard-periodic z. They are documented variants of Lars's setup.
- The ion-collision and hybrid runs use 259,776 fixed-position background
  macroparticles at n=1e16 m^-3, Ti=1 eV. The hybrid also has analytic Te=10 eV
  electron drag. Neither computes collective self fields.
- The kinetic runs use 259,776 mobile background ions and 259,776 mobile
  physical-mass electrons at ne=ni=1e10 m^-3, Ti=1 eV, Te=10 eV. They deposit charge,
  solve electrostatic fields, and execute six binary self/cross collision
  operators. Analytic electron stopping is absent in these cases.
- The profiler records binary collision execution in the collision runs and
  analytic stopping execution in the hybrid. This confirms dispatch, not an
  independently measured relaxation rate.

## Kinetic startup checks

Both kinetic runs start with exactly matched ion/electron spatial samples and
zero net background charge. Every background macroparticle remains in the
domain during the 10-ns test. Five full snapshots and every-step reduced
diagnostics are finite. Initial mean energies are approximately 1.500685 eV for
ions and 14.984536 eV for electrons, consistent with 3T/2 for the sampled Maxwellians.

- Base timestep: 5e-11 s, 200 steps. Half timestep: 2.5e-11 s, 400 steps.
- Combined initial Debye length: 0.0708799 m; cell size: 0.0625 m.
- Conservative electron gyro step:Omega_ce,max*dt=0.0812456.
- Final maximum |Ex|: 0.482505010 versus 0.482502934 V/m.
- Relative difference of this field observable: 4.30381e-06 (0.000430381 percent).
- Relative difference in final mean electron energy: 2.55434e-10.
- Minimum sampled electron distance from a z end: 0.221347 m.
- Remaining particle plus electrostatic-field energy changes by approximately
  1.57095e-05 relative to its initial value.

The last energy number excludes the prescribed static magnetic field and does
not explicitly subtract injected/escaped energy. The beam's physical weight is
tiny and no background losses occur here. It is an observed short-time energy
balance, not exact conservation or a complete open-system energy audit.

Halving dt doubles the number of injected macroparticle samples because the
source injects per timestep; their weights halve. Compare physical weights and
ensemble observables, not raw beam particle counts or individual trajectories.
This is a startup sensitivity check, not a mesh/particle-number convergence study.

The Ex response mainly comes from the initialized thermal column's evolution
and finite sampling. The inherited beam supplies only about 4000 physical ions/s
over the source area and represents only about 4e-5 physical ions over 10 ns.
Consequently this source does not meaningfully drive the background at its
present strength. The slight difference between the ion-only and hybrid beam
energies is not statistically established evidence of electron slowing.

## Geometry diagnosis and unresolved boundary

Two initial STL attempts are retained under `kinetic_initial20` and
`kinetic_grounded20` for debugging. The first hit a Poisson tolerance problem
with an almost cancelled RHS. Setting an absolute residual tolerance of 1e-10
V/m^2 handles that roundoff-level initial state. The second then stalled on a
nonzero RHS at a relative residual near 2.79e-3 with only one multigrid level.

The successful analytic-bore run uses the same strict 1e-8 relative tolerance and
1e-10 absolute tolerance, with four multigrid levels. No relaxed residual was
accepted as a successful solution. This implicates the STL/EB discretization
path, but does not identify a particular bad facet. The smooth cylindrical
inner wall approximates the original faceted bore and removes exterior fluid.

The installed source is upstream WarpX bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc.
The actual binary SHA256 is
`9cdadd080e1fa665c1c9a2f70a7e31ace1973f29d3943985d9a84d9a19d88085`.
It has no Periodic_ReflectParVel implementation. The located mirror field also
has transverse discontinuities at a standard periodic z join. Sampled particles
remain well away from that join in these runs; crossing physics is not tested.
See BOUNDARY_STATUS.md for the exact patch/branch information still needed.

## Files for the next session

For the kinetic configuration, open
`runs/kinetic_smoke200/inputs.txt`, followed by
`runs/kinetic_smoke200/warpx_used_inputs`, `metadata.json`, and `summary.json`.
The full data are in `diags/diag/openpmd/*.h5`; compact time histories are in
`diags/reduced/*.txt`. The half-dt and reservoir runs use the same layout.

`RESULTS.json` contains the measured values and individual boolean checks.
`test_workflow.py` has 8 passing composition/unit/scale tests; these do not replace
the GPU runs. Run `python report.py` to re-analyze these named runs and refresh
this report. Local runtime outputs are ignored by Git. Original input, STL,
and field files are preserved. No commits, pushes, or source rebuilds were made.

The next physics decision is the intended beam flux, density, and duration.
Increasing the density to the legacy 1e16 m^-3 while retaining this explicit
kinetic grid would be underresolved. Longer periodic trajectories require the
custom boundary/field mapping. These decisions affect what the experiment can
say, even though the electron/ion implementation already executes.
