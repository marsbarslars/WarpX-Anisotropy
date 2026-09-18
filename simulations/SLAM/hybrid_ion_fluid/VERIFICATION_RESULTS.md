# Electron-energy operator verification — 2026-09-18

The zero-drift thermal test closely follows the expected temperature-relaxation
curve. **The counterstream test exposes an energy-accounting problem in the
tested electron-ion relaxation path:** ion drift energy decreases without a
corresponding increase in electron internal energy. The zero-rate control does
not show that loss. Do not treat this operator as validated beam-to-electron
heating merely because the thermal test agrees with its reference.

These are small **numerical verification boxes, not a replacement mirror,
SLAM result, physical collision-rate prediction, or steady-state result**.
The physical mirror configuration and assets were not changed by these tests.

## What was tested

All three cases use an unmagnetized, fully periodic 0.5 m cube, 8×8×8 cells
(cell width 0.0625 m), collocated hybrid PIC, mobile proton particles, an
inertialess electron fluid with its energy equation enabled, and gamma=5/3.
The total proton/electron density is 10^16 m^-3. There are 4×4×4 markers per
cell **per ion species**: 32,768 markers in the thermal case; 65,536 across the
two counterstreams. Each species has its own Gaussian thermal velocities.

- **Thermal:** one proton species, nominal Ti=1 eV, Te=10 eV, no imposed drift.
- **Counterstream:** two proton species of density 5×10^15 m^-3 each,
  nominal Ti=Te=1 eV, drifting at +100 km/s and −100 km/s along x.
- **Counterstream control:** identical nominal initialization, but the
  electron-ion relaxation rate is zero.

The active tests impose a **constant, artificial nu=10^5 s^-1**. This value
isolates the numerical operator; it was not calculated for the laboratory
plasma. No ion-ion collisions, injection, wall, loss boundary, Joule heating,
or external B field is present. There are 100 steps of 1 ns, giving 100 ns,
and 10 field substeps per ion timestep. Runtime was about 10–11 seconds per
case on the existing GPU executable.

## Measured results

Energy comparisons start at **step 1, t=1 ns**, not at the pre-hybrid
initialization dump at step 0. All cases reach step 100 and retain their
particle counts. Temperature entries below are measured values, not the
nominal initialization parameters.

| Quantity | Thermal, nu=10^5/s | Counterstream, nu=10^5/s | Counterstream, nu=0 |
|---|---:|---:|---:|
| Electron T, 1 → 100 ns (eV) | 9.99820 → 9.82350 | 1.000000 → 0.999989 | 1.000000 → 1.000000 |
| Ion T, 1 → 100 ns (eV) | 1.00327 → 1.17759 | 1.001835 → 1.001851 | 1.001922 → 1.001922 |
| Maximum normalized thermal-reference misfit | **0.013237%** | Not applicable | Not applicable |
| Maximum absolute energy residual | **0.011107%** | **1.865099%** | **1.33089×10^-7%** |
| Signed final energy residual | −0.003518% | **−1.865099%** | **−1.33089×10^-7%** |
| Ion bulk-energy change (J) | −4.421×10^-9 | **−2.06112×10^-4** | −5.515×10^-13 |
| Electron internal-energy change (J) | −5.24797×10^-5 | **−3.16918×10^-9** | +2.636×10^-11 |

In the counterstream case, electron heating does not account for the
approximately 0.206 mJ removed from ion drift. Temperatures remain near 1 eV
while the opposing streams slow down. This is precisely why temperature
plots alone are insufficient to validate a slowing-down model.

### Definitions behind the numbers

For the uniform single-species, zero-drift thermal reference:

    Delta T(t) = Delta T(t_ref) exp[-4 nu (t - t_ref)]
    Te + Ti = constant

The plotted misfit is the measured minus reference temperature difference,
divided by the **initial temperature difference**, multiplied by 100. It is
not a percentage error in each individual temperature and is not a fitted
curve. The reference is anchored at the first initialized measured snapshot.

Ion temperature is obtained from physical-weighted particle velocities after
subtracting the **global mean of each species separately**. Opposing species
drifts are therefore not mistaken for thermal energy. This uniform-box
estimator is intentionally different from the mirror diagnostics' cell-local
pressure tensor. Native deposited species temperatures are also recorded in
the CSV; those use a local deposition estimator and need not equal the global
particle estimate at finite marker count.

The signed closed-box residual uses:

    E = sum_particles [w (gamma_particle - 1) m_p c^2]
        + integral [Pe / (gamma_e - 1)] dV
    residual = [E(t) - E(t_ref)] / E(t_ref)

Here gamma_particle is the relativistic particle factor, while gamma_e=5/3
is the electron heat-capacity parameter. Pe is in Pa. Native Te output is
Kelvin and is converted to eV. Electron temperatures in the plots are
density-weighted domain averages; electron energy is the volume integral of
pressure, not simply an unweighted mean temperature.

Magnetic energy is also evaluated and is zero in these three outputs. Electric
field energy is not included in this hybrid material-energy diagnostic;
this is not an electromagnetic PIC energy budget with displacement current.
No particle escape or source can explain the counterstream residual in this
closed periodic fixture.

## Interpretation and revision scope

The installed ion update is an Ornstein–Uhlenbeck drag/diffusion operation
toward the electron-fluid velocity:

    v_new = u_e + (v_old - u_e) exp(-nu dt) + thermal random increment

The electron-side relaxation sink uses the species thermal temperature
differences. In the inspected implementation, those thermal terms do not
provide the matching compensation for the relative bulk drift removed in
this counterstream test. This source reading and the paired numerical test
support the diagnosis for **this executable and configuration**; they do not
establish that every WarpX release, electron closure, or collision operator
has the same issue.

Verified runtime/source identity:

- WarpX revision: `bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc`.
- AMReX runtime: `26.08-18-g057940244648`.
- Executable SHA256:
  `9cdadd080e1fa665c1c9a2f70a7e31ace1973f29d3943985d9a84d9a19d88085`.
- Inspected source: `Source/FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.cpp`,
  electron sink near lines 1019–1110 and ion update near lines 1120–1300.
- No WarpX source patch or rebuild was performed for these tests.

The thermal result is encouraging **short-time agreement**, not exact
discrete conservation. Finite-marker stochastic variation, deposition,
interpolation, and operator splitting remain. Neither marker-count nor
timestep convergence was established by these three runs. A passing
zero-drift thermal check must not override the counterstream energy test.

**Practical decision:** keep this Qei/drag path experimental for moving-beam
heating. An evolving fluid electron-energy model can still be studied with
the unvalidated exchange path disabled and clearly labelled. Before enabling
it for quantitative beam heating, resolve the drift-energy accounting and
repeat thermal, counterstream, and timestep/marker-count verification.

## Files and reproduction

Existing results live below this directory:

- `runs/qei_thermal_100ns_20260918/`
- `runs/qei_counterstream_100ns_20260918/`
- `runs/qei_counterstream_control_100ns_20260918/`

Each contains `inputs.txt`, `manifest.json`, `warpx_used_inputs`, `warpx.log`,
`execution.json`, `verification_history.csv`, `verification_summary.json`, and
the slide-ready `energy_verification.png` / `energy_verification.svg`.
All three PNGs were visually inspected; labels, normalizations, and the
artificial-rate warning are visible and unobstructed. The panels autoscale:
the control's residual axis has a **10^-7 multiplier in percent**, whereas
the active-counterstream residual is of order 1 percent. Do not compare their
visual curve heights without reading those scales. Counterstream temperature
axes are also zoomed around 1 eV; their tiny variations are not strong heating.

From WSL Ubuntu, run the following. The `_repeat` directories must not already
exist; preparation deliberately refuses to overwrite a run. These commands
reproduce the parameters above and explicitly opt in to launching WarpX.

```bash
cd /mnt/c/Users/ryanv/Desktop/BOOGERAIDS/WIPPL/WarpX-Anisotropy/simulations/SLAM/hybrid_ion_fluid
PY=/home/ryanv/miniforge3/envs/warpx-gpu/bin/python

$PY verify_energy.py prepare runs/qei_thermal_100ns_repeat --case thermal --execute
$PY verify_energy.py analyze runs/qei_thermal_100ns_repeat

$PY verify_energy.py prepare runs/qei_counterstream_100ns_repeat --case counterstream --execute
$PY verify_energy.py analyze runs/qei_counterstream_100ns_repeat

$PY verify_energy.py prepare runs/qei_counterstream_control_100ns_repeat --case counterstream --nu 0 --execute
$PY verify_energy.py analyze runs/qei_counterstream_control_100ns_repeat

$PY -m unittest -v test_verify_energy
```

For CPU-only reanalysis of an existing run, use just `analyze` with its
existing directory. There are **10 passing unit tests** covering the thermal
reference, deck selections, parameter guards, weighted particle moments,
relativistic kinetic energy calculation, and signed residual definition.
