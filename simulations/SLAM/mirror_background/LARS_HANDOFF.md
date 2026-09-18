# Mirror-background temperature setup for discussion with Lars

Prepared for Ryan, 2026-09-14. This file is a local handoff draft, not a message
that has been sent. No Git push is performed by the simulation workflow.

## What we have

These are our **five workflow presets**, not five separate WarpX solvers. There
are four physical variants and one timestep sensitivity check.

| # | `--case` | Background and collision treatment | Self-consistent E | Can measure evolving electron temperature? |
| --- | --- | --- | --- | --- |
| 1 | `ion_background` | Fixed-position thermal protons; no collisions | No | No electron population |
| 2 | `ion_coulomb` | Fixed-position thermal protons; beam-background proton Coulomb collisions | No | No electron population |
| 3 | `hybrid` | Same ions/collisions plus analytic Maxwellian electron stopping | No | No: electron temperature is a prescribed bath parameter |
| 4 | `kinetic` | Mobile background protons and physical-mass electrons; all six self/cross-species Coulomb pairs | Yes, electrostatic | Yes; this is the temperature-plot case |
| 5 | `kinetic_half_dt` | Same physics as #4, with half dt and twice the steps | Yes, electrostatic | Yes; numerical check, not added physics |

Our name `hybrid` means the mixture of explicit ion collisions and an analytic
electron bath. It is **not** WarpX's separate hybrid-PIC/Ohm's-law field solver.
The analytic bath removes beam energy without evolving electron particles. Use
#4 or #5 to measure a dynamical electron temperature; do not add analytic electron
drag on top of the kinetic beam-electron collision pair.

The first three presets use `n=1e16 m^-3`; the kinetic presets use `n=1e10 m^-3`
to keep the electron/screening scale compatible with this explicit mesh. They
are not a controlled same-density comparison. Fixed-position ion particles can
still have their momenta changed by collisions; they are not a resetting thermostat.

## Current kinetic physical setup

- Three species: beam protons, background protons, and electrons with real masses.
- Initial background `Ti=1 eV`, `Te=10 eV`, `ni=ne=1e10 m^-3`.
- Background initial support: `r<0.75 m`, `0.25<z<4.75 m`; zero mean thermal drift.
- Prescribed mirror B from the located legacy HDF5 field; self-consistent
  electrostatic E from charge deposition and a Poisson solve. No evolved magnetic
  perturbations and no neutral transport/ionization.
- Pairwise Coulomb operators: beam-ion, beam-electron, electron-ion, ion-ion,
  electron-electron, beam-beam; pairwise Coulomb logarithm 10. No analytic drag.
- Pairwise energy/momentum correction is disabled, with highly unequal
  beam/background particle weights. Validate a full energy budget before
  interpreting heating; conservation is not enforced exactly for every pair.
- Domain: `-1<x<1`, `-1<y<1`, `0<z<5`, all in **metres**.
- Mesh `32x32x80`, cell width `0.0625 m`; dt `5e-11 s` (0.05 ns).
- Combined initial Debye length is about `0.0709 m`, only 1.13 grid cells:
  a basic scale check, not a grid-convergence demonstration.
- Grounded analytic cylindrical bore `r=0.9 m`, with absorbing particle wall;
  standard periodic z. This differs from the faceted STL/custom boundary input.

Lars's active continuous source is retained: an already-ionized proton injection
plane at `x=-0.2 m`, centered near `(y,z)=(0,1) m`, with `-0.1<y<0.1` and
`0.9<z<1.1 m`. It injects in the positive-x direction and also has positive-z
momentum. The source is **not at the midplane** and is not a 1000-particle pulse.
The inherited momentum spread produces approximately 9--16 eV protons.
The nominal flux `1e5 m^-2 s^-1` over `0.04 m^2` means about 4000 physical
protons/s; computational marker count is not physical beam strength.

## Exactly what "middle temperature" means

The geometric midplane is `z=2.5 m`, halfway along the 5-m box. We sample a
finite cylinder, not an infinitely thin plane and not the full cross-section:

```text
x^2 + y^2 < (0.25 m)^2
2.25 m < z < 2.75 m
```

Thus it is 0.5 m long and 0.5 m in diameter, centered on the axis. The output
filter does not constrain particle motion; particles can freely enter and leave.

For each species separately, the physical-weighted central temperature is
`T[eV] = m/(3*q_e) * <|v - <v>_w|^2>_w`. Mean motion is subtracted before computing
random velocity spread. This is a nonrelativistic second moment, not a Maxwellian
fit. Spatially varying flows inside the cylinder can still contribute spread.
Record timestamp, species temperature, component temperatures, mean flow, local
particle count and physical weight. No smoothing is used.

## How to change the model

Use `workflow.py --case ...` to generate a consistent input deck; WarpX itself
does not recognize our `--case` names. For example, in WSL from this folder:

```bash
/home/ryanv/miniforge3/envs/warpx-gpu/bin/python workflow.py prepare --case kinetic --boundary standard-periodic --steps 600 --temperature-every 1
```

`prepare` writes `runs/<name>/inputs.txt` and `metadata.json` without running.
Replace `prepare` with `run` to generate and execute a new run. For model #5,
replace `kinetic` with `kinetic_half_dt`; `--steps 600` then means 1200 actual
steps at half dt, covering the same 30 ns. For models #1--3, the same step count
does **not** imply the same duration because those presets have larger dt.

The generated text input is the actual WarpX configuration. You can run it
directly with `warpx.3d inputs.txt` after supplying the external field map path.
Do not switch a physical model by changing only its label: species, deposition,
particle pushing, collisions, density and dt must stay mutually consistent.
Our generator handles those grouped changes. Keep edited manual inputs as a
separate experiment with updated metadata rather than overwriting run provenance.

## Why not just extend to the beam-arrival time?

There are two different limits: duration and source strength. An axial straight
flight scale from z=1 to z=2.5 at about 3e4 m/s is 50 microseconds, not 50 ns.
That would require about one million electron-resolving steps at the present dt.
Actual mirror trajectories may not arrive on that simple estimate.

More time alone does not make the inherited source a useful heater. Even 50
microseconds of injection contributes only about 2.4 eV total, versus roughly
1.3e12 eV initial background thermal energy. To study heating, we need to choose
the source power/energy and background parameters deliberately, then distinguish
beam-caused changes from background evolution using a matched source-off case.

The current field map also does not join smoothly at the periodic z seam, and
Lars's custom boundary implementation is still unavailable locally. The initial
particle data give a maximum electron speed of about 7.43e6 m/s and a shortest
straight path-length / speed estimate of 36.5 ns to a 5-cm seam buffer. This is
an initial kinematic estimate, not a guarantee in the evolving electric field.
The unchanged-case extension is bounded at 30 ns, and particle/field diagnostics
are checked afterward. No arbitrary speed-changing boundary was introduced.

## What the existing 10-ns figure says

- Top-left: central ion temperature is almost constant near 1.007 eV; ion motion
  is much slower than electron motion. The offset from input 1 eV is finite sampling.
- Bottom-left: central electron temperature changes roughly 10.103 to 9.997 eV.
  That is a local second-moment transient, not demonstrated beam heating/cooling.
- Right panels: local sample counts. Electron count changes 3328 to 3204, while
  the global electron population stays 259776. Leaving the probe is not losing a
  particle from the device. Ion count remains 3328.
- Dashed lines are the input temperatures, not fitted equilibrium levels or
  thresholds. A curve crossing a dashed line is not proof of thermalization.
- No beam particles reach this central region during that 10-ns run.

## Completed extension: 30 ns

`runs/kinetic_temperature600_30ns_20260914` completed 600 steps in 123.61 s of
solver wall time, plus post-processing. There are 601 central samples, including
t=0. All integration checks and the 15 independent full-snapshot comparisons
passed. The latest 17 unit tests also passed, including the portable input check.

| Measurement | Initial | Final at 30 ns |
| --- | --- | --- |
| Central ion temperature | 1.006865 eV | 1.006864 eV |
| Central electron temperature | 10.102915 eV | 10.019087 eV |
| Central ion markers | 3328 | 3328 |
| Central electron markers | 3328 | 3221 |
| Beam markers in central region | 0 | 0 |

The electron temperature spans 9.957557--10.106345 eV and has a small bump around
23--24 ns rather than a sustained heating trend. Local electron counts span
3193--3328, but the global population remains 259776 for each background species.
Temperature changes are consistent with local population/early plasma evolution;
they do not isolate collisional beam heating without a control and energy budget.

The smallest electron distance to a z seam among the five full snapshots was
0.132789 m. This is a sampled clearance, not continuous tracking of every boundary
approach. The initial kinematic bound and saved snapshots support this bounded
startup, not an arbitrary extension to the beam-flight timescale.

No source flux, beam momentum, background density/temperature, wall, or boundary
was changed. Relative to the 10-ns input, only `max_step` and the ordinary full
snapshot interval changed. The central diagnostic remains every step.

Suggested message:

> We now have four physical presets plus a half-timestep kinetic check. The current
> run is the fully kinetic case: mobile background protons/electrons, Coulomb
> collisions, self-consistent electrostatic E, and prescribed mirror B. I extended
> the startup to 30 ns and recorded separate Ti/Te every step in r<0.25 m,
> 2.25<z<2.75 m. The temperature diagnostic passes its snapshot checks, but the
> inherited source is too weak and too far from this probe to establish beam
> heating here. I can send the standalone input and plot. Could we agree on beam
> power/location and the boundary treatment for a controlled heating comparison?
> The current variant uses an analytic cylinder and standard periodic z.

## What to send

Send this setup note, the chosen run's `inputs.txt`, `metadata.json`, temperature
PNG and compact summary. Share the HDF5 input field through agreed lab storage
with the checksum in `ASSETS.md`; update `particles.read_fields_from_path` on
Lars's machine. Large raw diagnostics stay out of Git. Source code and small
notes can be reviewed on the RyanZhu branch once pushed by Ryan's chosen workflow.

For convenient sharing, `inputs_kinetic_temperature_30ns.txt` is a standalone
copy of the executed 30-ns deck. Only the external field-map path is changed to
`example-femm-3d.h5` (and explanatory comments added). It does not depend on
`workflow.py` to run. Put the field beside it in a fresh output directory or
override that one path before launching the existing compatible 3D CUDA executable.
