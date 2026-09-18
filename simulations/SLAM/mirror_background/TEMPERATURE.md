# Central temperature: what is measured and why

## Physical question

Does the injected proton beam transfer directed energy into random motion of the
background ions and electrons? Separately track `Ti(t)` and `Te(t)` to investigate
where energy goes. Temperature also affects collisional rates, electron drag,
screening, and pressure, so it is an input to the physics as well as an output.

Do not call every change in local temperature beam heating. Particles can enter
and leave the measurement region, collective electric fields can redistribute
energy, and a finite computational sample has fluctuations. A heating study will
need an energy budget and a matched source-off comparison, not just a rising line.
No extra source-off campaign was run for this diagnostic addition.

## Definition and units

For each background species separately, over a finite cylinder centered on the
mirror axis (`r<0.25 m`, `2.25<z<2.75 m`):

```text
W = sum(w)
v_mean = sum(w*v) / W
T_j [eV] = m/q_e * sum(w*(v_j - v_mean_j)^2) / W
T [eV] = (T_x + T_y + T_z) / 3
```

`w` is the physical multiplicity represented by a computational macroparticle.
The denominator is physical weight, not marker count. The implementation uses
centered velocities for numerical stability. This is a distribution second
moment (no Bessel/sample-variance correction), not a Maxwellian fit. Spatially
varying flow within the cylinder can contribute to the unresolved velocity spread.

openPMD-viewer returns normalized momentum `u=p/(m*c)`. We first convert it with
`v=c*u/sqrt(1+|u|^2)`. The temperature definition itself is nonrelativistic, suitable
for these 1/10-eV backgrounds; the analyzer rejects particles faster than 0.05c.
The component temperatures are Cartesian, not parallel/perpendicular to B.

The distinction from kinetic energy is:

```text
mean kinetic energy [eV] = 3/2*T [eV] + m*|v_mean|^2/(2*q_e)
```

Thus a cold beam moving quickly has directed kinetic energy but almost no
temperature. Initial `Ti=1 eV` and `Te=10 eV` imply approximately 1.5 and 15 eV
of mean kinetic energy for zero-mean Maxwellian loading. A finite central sample
need not initially equal the requested temperature exactly.

Empty regions are recorded with missing temperature (JSON null / CSV blank),
not zero. Counts, physical weight, and effective count accompany the data.
Samples at adjacent times share particles and are not independent observations.

## Expectations for this startup case

The inherited source is already-ionized protons, not neutral-ionization physics.
Its nominal flux is `1e5 protons/(m^2 s)` across `0.04 m^2`, or about 4,000
physical protons/s. At roughly 12 eV this is only about `7.7e-15 W`. Computational
marker count is much larger than the injected physical weight and must not be
mistaken for the beam's physical strength.

In 10 ns, the source delivers only about `4.8e-4 eV` of total energy. The initial
kinetic background contains about `7.93e10` physical particles of each species
and roughly `1.3e12 eV` total thermal energy. Even depositing all that beam energy
in all background electrons would change their global temperature by only about
`4e-15 eV`. These are scale estimates, not a measured transfer rate.

The source is near `z=1 m`; the probe is at `z=2.5 m`. A simple straight axial
flight estimate using `v_z~3e4 m/s` is about 50 microseconds, compared with the
10-ns startup. Actual trajectories are magnetized. No direct beam-arrival heating
at the center should be inferred from this short run.

The analytic-electron `hybrid` case has a prescribed electron bath temperature;
it cannot produce a measured evolving electron temperature. The kinetic case
has actual mobile electron particles and can measure `Te`, but the present
duration and source strength do not establish beam-induced heating.

## Implementation and checks

The optional `center` full diagnostic uses WarpX's supported
`<diag>.<species>.plot_filter_function(t,x,y,z,ux,uy,uz)` and `fields_to_plot=none`.
It records all selected particles; no histogram approximation, random sample,
smoothing, or invented/interpolated timestamps are used. See the official
[diagnostic parameter documentation](https://warpx.readthedocs.io/en/latest/usage/parameters.html#running-cpp-parameters-diagnostics).
The installed checkout's `Docs/source/usage/parameters.rst` was also inspected.

`test_temperature.py` checks known temperature, bulk-drift invariance, unequal
weights, marker splitting, empty regions, normalized-momentum conversion,
nonrelativistic validity, and unchanged physical inputs when enabling the probe.
The analyzer independently selects the same region from the normal, unfiltered
full snapshots and checks agreement with the filtered diagnostic. These checks
validate measurement wiring; they do not demonstrate physical convergence.

## Verified run: 2026-09-14

`runs/kinetic_temperature200_20260914` completed 200 GPU steps at `dt=5e-11 s`,
covering 10 ns, in 47.10 s of solver wall time (post-processing is separate).
There are 201 central snapshots including the initial state. All 16 unit tests
passed, all integration checks passed, and 15 species/snapshot comparisons with
the ordinary unfiltered diagnostics agreed (three species at five times).

| Central population | Initial temperature | Final temperature | Initial/final marker count |
| --- | --- | --- | --- |
| Background protons | 1.006865 eV | 1.006865 eV | 3328 / 3328 |
| Electrons | 10.102915 eV | 9.997421 eV | 3328 / 3204 |
| Injected protons | Undefined (empty) | Undefined (empty) | 0 / 0 |

The initial temperatures differ slightly from input targets because this is a
finite local sample. Electrons move across the sampling boundary, changing which
particles are measured; their local temperature spans approximately 9.958--10.103
eV. This is not proof of whole-plasma cooling, beam heating, or equilibration.
The global electron marker count stayed at 259776: leaving the central probe is
not the same as leaving the device. No beam markers reached this probe.

The optional central raw diagnostics occupy about 123 MB for this test, and are
ignored by Git. Share the approximately 110-kB CSV, figure, and compact summary
when raw particle inspection is unnecessary. No long-duration heating run was made.

## Sharing with Lars

The unchanged kinetic setup has now also completed a 30-ns, 600-step extension:
`runs/kinetic_temperature600_30ns_20260914`. It produced 601 central samples and
passed the operational and independent-snapshot checks. Central Ti stays near
1.006864 eV; final central Te is 10.019087 eV, with no central beam markers.
See `LARS_HANDOFF.md` for the full model comparison, longer-run numbers, and a
standalone input file to share. Longer duration alone has not established heating.

Send the code/branch reference, the actual generated kinetic `inputs.txt`,
`metadata.json`, and the small temperature figure and summary. The generated
input contains a local absolute magnetic-field path: Lars must replace
`particles.read_fields_from_path` with his path to the agreed field map. Do not
send the parsed `warpx_used_inputs` dump as if it were the main handwritten deck.

The magnetic HDF5 map and bulky outputs belong in agreed lab storage, not this
Git commit. `ASSETS.md` records the field checksum and identity caveat. The kinetic
deck uses an analytic cylinder, so it does not require the original STL at runtime.

Report these limitations alongside the plots: analytic cylindrical wall instead
of the original faceted STL, explicit standard-periodic z instead of Lars's custom
boundary, the map's seam mismatch, and an intentionally short, low-density startup.
Ask Lars to agree on beam flux/energy, plasma density, and physical duration before
turning this into a heating study. Do not extend through the unresolved field seam
and present that as physical confinement.
