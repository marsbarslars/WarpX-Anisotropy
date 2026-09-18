# Electron-fluid implementation audit

Audit date: 2026-09-18. Scope: the installed WarpX source, not an assumption that online `latest` documentation matches the executable.

## Short explanation for Ryan, Lars, and Tony

We can let electron temperature evolve while keeping electrons as a fluid. The installed solver transports electron internal energy through an entropy variable and includes compression/expansion. Its pressure then acts back on the ions through Ohm's law. This is a meaningful improvement over fixing electron temperature at 10 eV, but it does not automatically include heat conduction, electron escape/sheaths, or beam-to-electron collisional heating.

The installed electron-ion relaxation option also damps ion motion toward the electron-fluid velocity. Source inspection found that the electron heat update accounts for the thermal temperature difference but does not explicitly receive the lost energy of an independently drifting ion population. We therefore keep this option **off in the production beam case**, rather than describe that run as energy-conserving electron-ion beam slowing. A separate, clearly labelled numerical verification case can test the option and show its residuals.

At the illustrative starting parameters, physical electron-ion thermal relaxation is slow: approximately 50 milliseconds for the electron-ion temperature difference to relax by a factor of e. A microsecond-scale run can show transport and ion dynamics without showing large collisional electron heating. Do not increase a physical collision coefficient merely to produce a more dramatic plot.

## Provenance and supported controls

Installed source root: `/home/ryanv/src/warpx` in WSL Ubuntu.

Audited clean source commit: `bd8c12ff8fc8d99682fef25cff92d99ea14f2bcc`.

Existing executable: `/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d`.

The references below are line numbers in that installed source. Source inspection establishes what that revision implements; a GPU run is still needed to establish executable/runtime behavior.

| Control or output | Meaning and exact installed-source evidence |
| --- | --- |
| `hybrid_pic_model.solve_electron_energy_equation = 1` | Enables the energy/entropy transport path. `Source/FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.cpp:87-100`, `:1319-1380`. |
| `hybrid_pic_model.gamma = 5/3` | Appropriate monatomic electron internal-energy factor for this energy-equation baseline. Do not retain the isothermal `gamma=1` when interpreting `Ue=Pe/(gamma-1)`. |
| `hybrid_pic_model.elec_temp = 10` | Initial/reference temperature is specified in eV and converted to joules at `HybridPICModel.cpp:131-132`; it is not a fixed-temperature constraint when the energy equation evolves it. |
| `hybrid_pic_model.n0_ref` | Required when gamma differs from 1, checked at `HybridPICModel.cpp:75-80`. |
| `hybrid_pic_model.include_joule_heating` | Optional resistive heating; off by default. Its coefficient must match the Ohm-law resistivity, not be invented independently. `HybridPICModel.cpp:100-107`, `:851-1009`. |
| `hybrid_pic_model.plasma_resistivity(rho,J,t)` | Global resistivity in ohm metres. This installed parser does **not** receive electron temperature as an argument. `HybridPICModel.cpp:82`; units in `Docs/source/usage/parameters.rst`. |
| `hybrid_pic_model.electron_ion_relaxation_rate(rho,Te,Ti,t)` | Optional rate in inverse seconds. `rho` is total positive ion charge density in C/m^3; `Te` and `Ti` are eV. Merely specifying it activates both the electron temperature sink and ion drag/diffusion. `HybridPICModel.cpp:114-125`. |
| `<ion species>.do_temperature_deposition` | Automatically enabled on every charged species when this relaxation option is requested: `Source/Particles/PhysicalParticleContainer.cpp:205-219`. Explicitly documenting it in a verification input is useful. |
| `hybrid_pic_model.joule_redirect_Te_threshold` | Artificial threshold redirection of electron Joule heating to ions. Leave off; it is not an electron-energy conservation correction. `HybridPICModel.cpp:109-113`. |
| Full-diagnostic `Te`, `Pe` | `Te` is **kelvin**, `Pe` is pascals. Convert `Te_eV = Te_K * kB / qe`. `HybridPICModel.cpp:175-194`; `Docs/source/usage/parameters.rst:4468`. |

The production transport upgrade should use the energy flag, gamma=5/3, and zero resistive/Joule/relaxation source terms unless those terms are separately selected and verified. Existing ion-ion Coulomb collisions are a separate operator and do not turn on electron-ion energy exchange.

### Fluid-velocity diagnostic registry workaround

The first runtime attempt to export `hybrid_electron_velocity_fpx` (and y/z) stopped during diagnostic initialization, before time advancement. Although `FullDiagnostics.cpp:949-954` appears to support a vector base name followed by a direction, its condition calls the **scalar** `has(base, level)` rather than a vector-component lookup. This electron-velocity field is allocated only as vector components, so that condition fails.

The input instead requests explicit component names as scalar registered-field references:

```text
diag.additional_fields_to_plot = "hybrid_electron_velocity_fp[dir=x]" "hybrid_electron_velocity_fp[dir=y]" "hybrid_electron_velocity_fp[dir=z]"
```

`Source/ablastr/fields/MultiFabRegister.cpp:665-686` constructs a vector component's complete key as `base[dir=x][level=0]`; scalar lookup appends only `[level=0]`. Therefore the explicit `[dir=x]` name reaches the real component through the generic scalar diagnostic fallback at `FullDiagnostics.cpp:1009-1010`, without changing WarpX source. The manifest records these raw output names in `electron_velocity_output_fields`. Read them as three scalar fields, in m/s, rather than as an electron-particle distribution. **GPU confirmation completed:** the 10-ns retry, 100-ns adaptive test, and both final 1-microsecond runs accepted these names; analysis successfully read all three velocity arrays. Unit tests alone were not used as evidence of runtime acceptance.

## What the transported fluid represents

The electron density is inferred from quasi-neutrality. The electron velocity is

`ue = (Ji - Jplasma) / rho`,

implemented at `HybridPICModel.cpp:709-723`. The solver transports an entropy variable proportional to

`Ke = Te * ne^(1-gamma)`

using numerical Lagrangian markers, then reconstructs temperature and `Pe = ne kB Te`. These markers are an advection device, **not kinetic electrons or samples of an electron velocity distribution**. See `HybridPICModel.cpp:731-848`, `:1289-1315`, `:1319-1380`.

The corresponding source-free fluid equation is

`dUe/dt + div(Ue ue) + Pe div(ue) = 0`, with `Ue = Pe/(gamma-1)`.

Thus temperature can change through transport and compression even with electron-ion collisions disabled. No conductive heat flux is included in this transport equation. Electron pressure is scalar/isotropic; ion parallel and perpendicular pressures remain distinct measurable particle moments.

The [official hybrid-model explanation](https://warpx.readthedocs.io/en/latest/theory/models_algorithms/kinetic_fluid_hybrid_model.html) describes this approximation and entropy-transport formulation. Its source-dependent collision details must be interpreted with the version discrepancy below.

## Why the installed relaxation option is not yet approved for beam energy accounting

### Electron update

For each ion population s, the electron sink is

`Te_new = Ti_s + (Te_old - Ti_s) exp[-3(gamma-1)(ns/ne) nu dt]` for singly charged ions.

This is the actual kernel at `HybridPICModel.cpp:1100-1110`. It depends on the **thermal** ion temperature, not its directed kinetic energy. Ion temperature deposition subtracts the species mean velocity; see `PhysicalParticleContainer.cpp:2018-2230`. That deposition also uses a sample-count correction `n/(n-1)`, while a particle energy budget uses population moments; the quantities are not exactly interchangeable at small marker counts.

### Ion update

The actual installed stochastic update is

`v_new = ue + (v_old-ue) exp(-nu dt) + sqrt[(kB Te/mi)(1-exp(-2 nu dt))] R`,

where R has three independent standard-normal components. See `HybridPICModel.cpp:1270-1283`. Therefore it damps an ion population's bulk velocity relative to **electrons**, not merely the population's own thermal spread.

No term returning `mi |ui_s-ue|^2` drift energy to the electron fluid appears in this update or the calling sequence at `HybridPICModel.cpp:1380-1510`. The Joule source instead uses `Jplasma/(e ne)` for every species and reduces to global `eta Jplasma^2`: `HybridPICModel.cpp:967-1002`. It cannot recover arbitrary counterstreaming populations' lost drift energy when their total current cancels.

### Decisive counterexample

Take two equal-density proton populations in a uniform periodic box, with equal thermal temperatures `Ti_1=Ti_2=Te`, opposite drifts `+U` and `-U`, zero magnetic field, and zero resistivity. Total ion current and electron velocity are zero.

The electron temperature sink is initially zero. In expectation, the ion update preserves each population's thermal variance but damps its mean drift by `exp(-nu dt)`. The lost directed energy density is

`sum_s (ns mi Us^2/2) [1-exp(-2 nu dt)]`.

There is no matching electron gain in the installed algorithm. Reducing dt does not eliminate this underlying missing transfer. This is a source-derived counterexample, not a claim that a production GPU run has already measured its magnitude.

### Even zero-drift thermal exchange is not exactly conservative per realization

The electron step is done first; the ion stochastic step then uses the updated electron temperature. In the uniform equal-density, gamma=5/3, constant-rate, infinite-marker idealization, let `a=exp(-2 nu dt)` and `D=Te-Ti`. The sequential updates give

`(Te_new+Ti_new)-(Te+Ti) = -D (1-a)^2`.

This splitting error converges as dt decreases. Finite-marker random kicks add sampling fluctuations; interpolation and temperature-moment definitions add further differences. The installed documentation's phrase “conserves energy exactly” should therefore not be used as an unqualified description of this executable.

## Physically scaled hydrogen rate and units

For nonrelativistic, weakly coupled hydrogen and a Maxwellian electron bath, with ion drift and ion thermal speed small compared with electron thermal speed, the small-speed ion friction coefficient is

`nu_ie = [4 sqrt(2 pi)/3] ne e^4 lnLambda sqrt(me) / [(4 pi epsilon0)^2 mp (e Te_eV)^(3/2)]`

or numerically

`nu_ie [s^-1] = 1.58283e-15 * ne [m^-3] * lnLambda / Te_eV^(3/2)`.

This follows from electron-ion momentum relaxation and momentum conservation: `nu_ie=(me/mp) nu_ei_momentum` for equal electron/proton density. The [NRL Plasma Formulary](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Formulary_2019.pdf) lists `nu_e = 2.91e-6 ne lnLambda Te^(-3/2)` with density in **cm^-3**; converting density and multiplying by the electron/proton mass ratio gives the rounded coefficient above. The unrounded number here is evaluated from SI constants. [PlasmaPy's Maxwellian-collision documentation](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.collisions.frequencies.MaxwellianCollisionFrequencies.html) states the slow-relative-flow and Maxwellian applicability conditions.

For hydrogen, the installed rate parser's OU/thermal conventions require **nu_ie**, not the roughly 1836-times-larger electron momentum frequency. With gamma=5/3 and a single equal-density proton species, `dTi/dt = 2 nu_ie (Te-Ti)` and `dTe/dt = -2 nu_ie (Te-Ti)` in the continuum zero-drift limit; the temperature difference decays at `4 nu_ie`.

Illustrative starting parameters: `ne=1e16 m^-3`, `Te=10 eV`, `Ti=1 eV`, **assumed** `lnLambda=10` (a model choice, not a fitted measurement):

| Quantity | Value |
| --- | --- |
| Ion friction coefficient nu_ie | 5.0053 s^-1 |
| Electron momentum-relaxation frequency | 9190.6 s^-1 |
| Temperature-difference e-folding time, 1/(4 nu_ie) | 49.95 ms |
| Initial background-ion warming slope | 90.10 eV/s |
| Linearized ion temperature rise over 200 microseconds | 0.0180 eV |
| Linearized ion temperature rise over 50 ns | 4.50e-6 eV |
| Speed of a 12 eV proton | 4.79e4 m/s |
| Electron most-probable thermal speed at 10 eV | 1.88e6 m/s |

The speed ratio is about 0.026, consistent with the slow-ion approximation. Here electrons start hotter than background ions: the initial thermal exchange would **cool electrons and warm background ions**. It must not be presented automatically as beam heating of electrons.

A parser expression illustrating these conventions is:

```text
# NOT enabled in the production beam run; source-energy caveat above applies.
# rho/qe is electron density for a quasi-neutral, singly charged hydrogen plasma.
# hybrid_pic_model.electron_ion_relaxation_rate(rho,Te,Ti,t) = "1.58282966e-15*(rho/1.602176634e-19)*10/max(Te,1.e-3)^1.5"
```

For appreciable ion temperatures relative to `(mp/me) Te`, large beam speeds relative to electron thermal speed, other ion masses/charges, non-Maxwellian electrons, or strong coupling, this approximation and coefficient must be revised. The numerical Te floor is not a physical cold-electron model.

Using the same initial collision scaling gives a resistivity scale `me nu_ei_momentum/(ne e^2) = 3.26e-5 ohm m`. A single constant value is only a frozen-temperature estimate: it is not a temperature-responsive resistivity for an evolving fluid. Moreover, the installed global resistivity parser lacks Te as an argument. Leave eta and Joule heating zero in the transport baseline rather than quietly call a constant a fully evolving Spitzer model.

## Installed source versus current online documentation

As inspected on 2026-09-18, the online hybrid theory page says that relaxation kicks occur about **each ion population's bulk velocity**, preserving bulk momentum. It also documents a per-species resistivity overlay with a Te-dependent parser.

The installed clean revision instead kicks about the **electron velocity**, as shown above, and has only the global `plasma_resistivity(rho,J,t)` parser in this model. No per-species resistivity implementation was found in its hybrid source. Do not copy those newer online parameters into this input and assume they are active. This discrepancy is an important reason to record the source commit and executable provenance with each result.

## Boundaries and remaining validity limits

No hard 3D Cartesian embedded-boundary prohibition was found; the constructor explicitly rejects the energy equation only in RCYLINDER/RSPHERE geometries at `HybridPICModel.cpp:95-100`. This is not proof that every vessel boundary is physically represented.

`Source/Fluids/QdsmcParticleContainer.cpp:122-129` creates one entropy marker for every cell. This file has no embedded-boundary intersection or cut-cell-volume treatment. Its `PushX` at `:349-365` clamps entropy markers at nonperiodic box boundaries rather than modelling escaping-electron energy; periodic directions wrap. Below the density floor, electron velocity is set to zero (`HybridPICModel.cpp:709`), and the entropy conversion uses a floored density to provide an insulating halo (`:766-773`). Thus:

- Interior, whole-cell central and near-source diagnostics are the cleanest first measurements.
- Do not interpret the floor/halo or marker clamp as physical electron wall loss, sheath physics, or a conductive heat-flux boundary.
- Audit cut-cell contributions before claiming a complete vessel energy budget.
- Check electron **fluid advection** distance `max(|ue|) dt/dx` as well as field substep stability; fluid electrons remove kinetic electron orbits, not every fast numerical timescale.
- The model neglects electron kinetic tails and heat conduction. Small electron mass motivates the reduction, but does not establish that these omitted effects are negligible for future SLAM confinement/steady-state questions.

## Verification and honest misfit plots

1. **Transport:** use a controlled periodic compression with initially uniform entropy. Compare measured `Te` to `Te0 (ne/n0)^(gamma-1)`, and report density-weighted residuals over valid populated cells. If initial entropy is not uniform, this pointwise formula is not a valid universal reference; compare advected entropy instead.
2. **Internal diagnostic identity:** check `Pe = ne kB Te` and the kelvin-to-eV conversion. This validates field reading/units; it is **not independent verification of the energy equation** because Pe is constructed from those quantities.
3. **Thermal-exchange verification only:** uniform periodic box, no beam drift, no fields/resistivity, single proton species, chosen constant nu. Compare `Te-Ti` against `exp(-4 nu t)` and plot total electron internal plus ion kinetic energy residual. Repeat at half dt and/or higher marker count. An artificially large benchmark nu must be labelled numerical, not the laboratory collision rate.
4. **Directed-energy regression:** use the equal counterstreaming populations described above; report drift damping and total energy residual. A zero-drift relaxation benchmark passing does not settle this issue.
5. **Production mirror:** plot evolving ion/electron temperatures, pressure, density, separate beam/background distributions, and timestep-comparison residuals. Do not claim electron-ion collisional heating, steady state, or accurate wall heat loss while those terms remain disabled or unvalidated.

The installed official relaxation benchmark lives at `Examples/Tests/ohm_solver_electron_energy_eq/inputs_test_2d_ohm_solver_electron_energy_picmi.py:493-565`, with analysis in `analysis_qei.py`. It is a zero-drift verification case; its default tolerances allow a 5% rate error and 2% thermal-energy drift. It does not test the independent-beam drift-energy issue.
