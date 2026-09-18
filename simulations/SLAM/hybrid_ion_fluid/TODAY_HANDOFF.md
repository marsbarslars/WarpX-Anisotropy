# Today's mirror-fluid work — 18 September 2026

## Start here

The new setting is **kinetic protons + an evolving, massless electron fluid**.
Electron temperature is no longer prescribed to stay at 10 eV: internal energy
is transported and responds to compression/expansion. Electron pressure feeds
back on ion motion through the hybrid field solver.

This is still the existing **analytic periodic mirror test**, with the inherited
old mirror source. It is not the original FEMM field map and not the SLAM vessel.
No new physical geometry was introduced today; no vessel integration, rebuild,
merge, commit, or push was performed.

**Important limit:** electron-ion collisional beam slowing/heating is OFF.
Testing the installed relaxation operator exposed loss of independent ion drift
energy without a matching electron gain. See `VERIFICATION_RESULTS.md` for the
measured counterstream/control test, and `ELECTRON_MODEL_AUDIT.md` for equations
and installed-version source evidence. This is not a claim about every WarpX
version. Ion-ion collisions also remain off in today's transport baseline;
the previously completed 50-ns ion-collision test is preserved.

## Verified completion

Both final CUDA runs reached **1 microsecond** with exit status zero and finite
diagnostics. The 5-ns/200-step case took 126.76 s; the 2.5-ns/400-step reference
took 134.42 s. Each has 11 full snapshots, including initialization. All **63
hybrid/model/plot unit tests** and **9 existing temperature tests** passed.

Final measurements from the 5-ns run:

| Quantity | Center | Near source |
|---|---:|---:|
| Volume-mean electron Te | 9.999428 eV | 9.976211 eV |
| Electron cell Te range | 9.91672–10.09268 eV | 9.87368–10.08003 eV |
| Ion T_parallel | 0.962001 eV | 0.904352 eV |
| Ion T_perp | 0.950129 eV | 0.960835 eV |
| Ion P_parallel | 1.540197 mPa | 1.447255 mPa |
| Ion P_perp | 1.521190 mPa | 1.537646 mPa |
| Fluid Pe | 16.020188 mPa | 15.928436 mPa |
| Ion density from cell counts | 9.992877e15/m³ | 9.988426e15/m³ |
| Fluid density from deposited charge | 9.999231e15/m³ | 9.964562e15/m³ |
| Background markers in region | 11,224 | 1,726 |
| Beam markers in region | 0 | 4,090 |

All 974,400 background markers remain in the whole domain; particles can move
between diagnostic cells without being lost. The domain contains 4,108 beam
markers representing only **0.00401172 physical protons**. No beam reaches the
central diagnostic region. Initial measured ion scalar temperatures were
0.958465 eV (center) and 0.944566 eV (near source), despite nominal 1-eV loading:
finite sampling and per-cell flow subtraction account for that initial offset.

The maximum matched-sample timestep difference across background-ion T/P/n is
**0.033067% of the documented initial reference scale**. For fluid Te/Pe/ne it is
**0.000464594%**. These are small numerical differences, not accuracy against
experiment or relative errors in the much smaller heating signal. The largest
ion discrepancy is center P_parallel at 0.7 microseconds; hard cell selections
can jump when a particle crosses a selected-cell boundary. Two timesteps do
not establish a convergence order, grid convergence, or marker convergence.
The beam injection uses a fixed marker count per source cell per step, so
changing dt also changes beam sampling/weight; beam PDFs are not included in
this background-moment convergence comparison.

The saved electron-advection CFL proxy peaks at 0.04024 (5 ns) and 0.01886
(2.5 ns). It is a cell-averaged proxy, not a strict internal CFL bound. Neither
diagnostic region intersects the embedded boundary or drops below the density
floor. The global final electric-field maximum is about 550 V/m, concentrated
in a non-equilibrium startup configuration; successful interior diagnostics do
not validate wall/edge physics.

Six slide-ready PNG/SVG figure pairs were generated and visually checked. The
temperature curves use clearly labelled separate ion/electron scales. Do not
call the nearly flat curves steady state: the horizon is short and the source
is intentionally tiny. Today establishes a usable transport/diagnostic baseline,
not the final physical heating model.

## Input reading order

Open `runs/energy_transport_1us_20260918/inputs.txt`: this is the actual executed
WarpX input, not merely a planning template. Read it in this order:

1. **Domain and mesh:** x,y = -1 to 1 m; z = 0 to 5 m; 32×32×80 cells.
   Each cell is 0.0625 m wide. One level, no adaptive mesh refinement.
2. **Solver and electron model:** `algo.maxwell_solver = hybrid`, gamma=5/3,
   `solve_electron_energy_equation = 1`. Initial Te=10 eV, reference n=1e16/m³,
   density floor=1e15/m³. No electron particle species is created.
3. **Time controls:** particle/fluid step 5 ns, 200 steps = 1 microsecond.
   RKF45 adapts magnetic-field substeps, not the mesh or main particle step.
   Full snapshots every 20 steps = 100 ns, plus initialization/final output.
4. **Applied mirror field:** B=curl(A), smooth periodic z dependence.
   On-axis Bmin=0.001487675 T at z=2.5 m; Bmax/Bmin=3.10951921.
5. **Walls and wrap:** analytic radius-0.9-m bore; absorbing ion walls and
   ordinary periodic z boundaries. A wrap is not a loss or Lars's custom boundary.
6. **Source:** continuous already-ionized protons at x=-0.2 m, directed +x,
   y=-0.1..0.1 m, z=0.9..1.1 m. Normalized initial momenta ux,uz=0.0001..0.00011,
   uy=-0.000095..0.000095, corresponding roughly to 9.38–15.58 eV, about 12 eV mean.
   Input flux=1e5/m²/s over 0.04 m², or 4000 represented physical protons/s.
   This is an extremely weak test source, not a useful laboratory heating power.
7. **Background:** mobile, depositing kinetic protons, initial Gaussian Ti=1 eV,
   n=1e16/m³ inside r<0.75 m over the full periodic z extent, 27 markers/cell.
   This sharp-edged initial column is not a force-balanced equilibrium.
8. **Diagnostics:** both ion species plus E, B, charge density, fluid Te/Pe and
   electron-fluid velocity. `Te` in the raw output is kelvin; analysis converts
   to eV. The first pressure dump is uninitialized and excluded from fluid curves.

## What the diagnostics measure

The two fixed, whole-cell regions are:

| Region | Cells and physical extent | Purpose |
|---|---|---|
| Center | 416 cells; cell centers r<0.25 m; z=2.25..2.75 m, k=36:44 | Interior low-field plasma, separate from injection |
| Near source | 64 cells; i=12:16,j=14:18,k=14:18; x=-0.25..0 m,y=-0.125..0.125 m,z=0.875..1.125 m | Source and immediate downstream region |

Indices are zero-based and upper bounds exclusive. `sampling_cells.json` saves
every cell and the actual summed cell volume. These masks are not transferable
to Lars's new mesh without redesign.

- Ion velocity distributions use local B and physical particle weights, with
  beam and background kept separate. A fluid has no measured electron velocity PDF.
- Ion T_parallel/T_perp and P_parallel/P_perp use velocities relative to each
  cell's own mean flow. Thus directed beam motion is not mislabeled temperature.
- Fluid Te is shown as both a volume mean and <Pe>/(e<ne>); spatial averaging
  means these need not be identical. Electron pressure is scalar in this model.
- Density, computational marker counts, represented physical-particle counts,
  histogram captured weight, and sparse-cell flags remain available in the data.
- Misfit means **difference from the half-timestep run**, scaled by stated
  initial quantities. It is not a fit to experiment or proof of physical accuracy.

## Plots and data locations

Final run: `runs/energy_transport_1us_20260918/`.
Reference: `runs/energy_transport_halfdt_1us_20260918/`.
Slide-ready files: final run's `presentation/`, PNG for easy insertion and SVG
for editable vector figures. No slide deck or new animation is generated.

- `01_temperature`: ion parallel/perpendicular T and electron-fluid Te.
- `02_pressure_density`: ion/fluid pressure and background density.
- `03_particle_inventory`: marker counts versus physical population.
- `04_numerical_misfit`: background-ion timestep comparison.
- `05_velocity_space`: background at the center and beam near injection.
- `06_fluid_numerical_misfit`: electron-fluid timestep comparison.

Underlying data are `ion_history.csv`, `electron_fluid_history.csv`,
`velocity_distributions.npz`, and `diagnostic_audit.json`. Actual model inputs,
binary/source-input hashes, runtime log and exit status are saved with each run.
Open `diagnostics.py` to inspect the moment definitions; `setup_case.py` generates
the input, and `presentation_plots.py` makes the figures.

## What to say to Lars and Tony

We have moved from fixed-temperature electrons to an evolving electron-fluid
energy model, keeping kinetic background and injected protons. Defined central
and near-source cell regions now report velocity distributions, anisotropic ion
temperature/pressure, fluid temperature/pressure, density and population counts.
The plots are startup/diagnostic results, not steady-state heating results.

The separate energy-operator audit is worth discussing before adding electron
drag: its thermal relaxation test agrees closely, but a counterstream test in
our installed revision loses directed ion energy without heating the electrons
by the corresponding amount. We therefore have not enabled that option in the
mirror or claimed a complete slowing-down model.

For a physical heating study, agree next on beam flux/energy, electron collision
and heat-transport closure, source/loss/energy balance, and a beam-off control.
Those are follow-up decisions, not extra implementation work scheduled today.

## Vessel and sharing status

Lars's fetched `origin/main` commit `50d81ac` contains matching updated STL
objects for its Mirror and SLAM directories (blob `a7694e6a...`, 1,237,084 bytes).
The local working SLAM STL is still the older 387,384-byte version: the updated
asset is available in fetched Git history, not silently installed here. It is
not used for today's tests.

Keep code, tests, and reviewed Markdown notes in Git. Runs/HDF5 dumps and generated
plots stay ignored; share selected figures separately or deliberately choose a
small reviewed set later. No raw output or executable was added to Git today.

Stop here after these verified results; the next activity is reading the
input together, not another simulation campaign.
