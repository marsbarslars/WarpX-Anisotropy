#!/usr/bin/env python3
"""Generate the hybrid collisional-mirror smoke test and production matrix.

The model deliberately stays test-particle-like:

* beam protons move in prescribed B;
* thermal background protons keep fixed spatial sampling but participate in
  pairwise momentum exchange;
* both species skip charge/current deposition, so no collective E response;
* beam/background-ion pitch scattering uses pairwise Coulomb collisions;
* electron drag uses WarpX's fixed analytic background-stopping model.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PARAMETERS = json.loads((ROOT / "parameters.json").read_text(encoding="utf-8"))
RUNS = ROOT / "runs"

Q_E = 1.602176634e-19
M_P = 1.67262192595e-27
C = 299792458.0


def density_tag(density: float) -> str:
    exponent = int(round(math.log10(density)))
    if not math.isclose(density, 10.0**exponent, rel_tol=1.0e-12):
        raise ValueError("Density tags currently require exact powers of ten")
    return f"1e{exponent}"


def ion_thermal_u_std(temperature_eV: float) -> float:
    """One-component nonrelativistic thermal std, sqrt(kT/m)/c."""
    return math.sqrt(temperature_eV * Q_E / M_P) / C


def input_text(*, angle_deg: int, density: float, max_step: int, dt_s: float,
               collision_model: str, beam_diag_interval: int,
               background_diag_interval: int, beam_particles: int,
               random_seed: int) -> str:
    domain = PARAMETERS["domain"]
    birth = PARAMETERS["beam_birth"]
    angle_rad = math.radians(angle_deg)
    u0 = PARAMETERS["beam_u0"]
    ion_u_std = ion_thermal_u_std(PARAMETERS["background_ion_temperature_eV"])
    ppc = PARAMETERS["background_particles_per_cell_each_dim"]
    cells = domain["cells"]
    center = birth["center_m"]
    rms = birth["rms_m"]

    collision_names = []
    collision_sections = []
    if collision_model in {"hybrid", "ion_only"}:
        collision_names.append("beam_ion_coulomb")
        collision_sections.append(f"""beam_ion_coulomb.type = pairwisecoulomb
beam_ion_coulomb.species = beam_protons background_ions
beam_ion_coulomb.CoulombLog = {PARAMETERS['coulomb_log']}
beam_ion_coulomb.ndt_supercycle = 1""")
    if collision_model in {"hybrid", "electron_only"}:
        collision_names.append("beam_electron_stopping")
        collision_sections.append("""beam_electron_stopping.type = background_stopping
beam_electron_stopping.species = beam_protons
beam_electron_stopping.background_type = electrons
beam_electron_stopping.background_mass = m_e
beam_electron_stopping.background_density = n_bg
beam_electron_stopping.background_temperature = Te_eV*q_e/kb""")
    if collision_model == "none":
        collision_block = "# Collisions disabled for the matched mirror-only control."
    elif collision_names:
        collision_block = (
            f"collisions.collision_names = {' '.join(collision_names)}\n\n"
            + "\n\n".join(collision_sections)
        )
    else:
        raise ValueError(f"Unsupported collision model: {collision_model}")
    electron_constant = (
        f"my_constants.Te_eV = {PARAMETERS['background_electron_temperature_eV']:.17g}"
        if collision_model in {"hybrid", "electron_only"}
        else ""
    )

    return f"""##########################
# HYBRID COLLISIONAL MIRROR
##########################
# Tutorial-scale beam: {PARAMETERS['beam_energy_eV']:.6f} eV, born at B-min.
# Background ions retain a kinetic velocity distribution and scatter the beam,
# but their positions are fixed so an absorbing boundary cannot drain the bath.
# Electrons are an analytic stopping bath and are not kinetic particles.
# Charge/current deposition is disabled for both kinetic species.

my_constants.Lx = {domain['Lx_m']}
my_constants.Ly = {domain['Ly_m']}
my_constants.Lz = {domain['Lz_m']}
my_constants.dt = {dt_s:.17g}
my_constants.Nbeam = {beam_particles}
my_constants.u0 = {u0:.17g}
my_constants.alpha_rad = {angle_rad:.17g}
my_constants.n_bg = {density:.17g}
{electron_constant}
my_constants.ui_std = {ion_u_std:.17g}

############
# NUMERICS #
############
geometry.dims = 3
geometry.prob_lo = -0.5*Lx -0.5*Ly 0.0
geometry.prob_hi =  0.5*Lx  0.5*Ly Lz
amr.n_cell = {cells[0]} {cells[1]} {cells[2]}
amr.max_level = 0
max_step = {max_step}
warpx.const_dt = dt

algo.particle_shape = 1
warpx.do_electrostatic = labframe
warpx.grid_type = collocated
warpx.serialize_initial_conditions = 1
warpx.use_filter = 0
warpx.random_seed = {random_seed}

##############
# BOUNDARIES #
##############
boundary.field_lo = pec pec pec
boundary.field_hi = pec pec pec
boundary.particle_lo = absorbing absorbing absorbing
boundary.particle_hi = absorbing absorbing absorbing

#############
# PARTICLES #
#############
particles.species_names = beam_protons background_ions

# Minority beam: identical spatial seed and birth geometry to the collisionless scan.
beam_protons.species_type = proton
beam_protons.do_not_deposit = 1
beam_protons.initialize_self_fields = 0
beam_protons.injection_style = gaussian_beam
beam_protons.x_m = {center[0]}
beam_protons.y_m = {center[1]}
beam_protons.z_m = {center[2]}
beam_protons.x_rms = {rms[0]}
beam_protons.y_rms = {rms[1]}
beam_protons.z_rms = {rms[2]}
beam_protons.npart = Nbeam
beam_protons.q_tot = q_e*Nbeam
beam_protons.momentum_distribution_type = constant
beam_protons.ux = u0*sin(alpha_rad)
beam_protons.uy = 0.0
beam_protons.uz = u0*cos(alpha_rad)

# Fixed-position Maxwellian proton reservoir. Momentum is still exchanged by collisions.
background_ions.species_type = proton
background_ions.do_not_deposit = 1
background_ions.do_not_push = 1
background_ions.initialize_self_fields = 0
background_ions.injection_style = NUniformPerCell
background_ions.num_particles_per_cell_each_dim = {ppc[0]} {ppc[1]} {ppc[2]}
background_ions.profile = constant
background_ions.density = n_bg
background_ions.momentum_distribution_type = maxwellian
background_ions.maxwellian_u_std_distribution_type = constant
background_ions.ux_std = ui_std
background_ions.uy_std = ui_std
background_ions.uz_std = ui_std
background_ions.maxwellian_u_mean_distribution_type = constant
background_ions.ux_mean = 0.0
background_ions.uy_mean = 0.0
background_ions.uz_mean = 0.0

##############
# COLLISIONS #
##############
{collision_block}

##########
# FIELDS #
##########
particles.B_ext_particle_init_style = read_from_file
particles.read_fields_from_path = ../../example-femm-3d.h5

###############
# DIAGNOSTICS #
###############
diagnostics.diags_names = beam_diag background_sample

beam_diag.diag_type = Full
beam_diag.format = openpmd
beam_diag.fields_to_plot = none
beam_diag.intervals = {beam_diag_interval}
beam_diag.species = beam_protons
beam_diag.beam_protons.variables = x y z ux uy uz w

background_sample.diag_type = Full
background_sample.format = openpmd
background_sample.fields_to_plot = none
background_sample.intervals = {background_diag_interval}
background_sample.species = background_ions
background_sample.background_ions.variables = x y z ux uy uz w
background_sample.background_ions.random_fraction = {PARAMETERS['background_diag_random_fraction']}

warpx.reduced_diags_names = particle_energy particle_momentum particle_number
particle_energy.type = ParticleEnergy
particle_energy.intervals = {beam_diag_interval}
particle_energy.precision = 18
particle_momentum.type = ParticleMomentum
particle_momentum.intervals = {beam_diag_interval}
particle_momentum.precision = 18
particle_number.type = ParticleNumber
particle_number.intervals = {beam_diag_interval}
particle_number.precision = 18
"""


def run_script(case_name: str) -> str:
    return f"""#!/usr/bin/env bash
set -u

WARPX=/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d
if [[ ! -x \"$WARPX\" ]]; then
    echo \"WarpX executable missing: $WARPX\" >&2
    exit 2
fi

echo \"Running hybrid collisional mirror case: {case_name}\"
\"$WARPX\" inputs 2>&1 | tee run.log
status=${{PIPESTATUS[0]}}
echo \"WarpX exit code: $status\"
exit \"$status\"
"""


def write_case(case_name: str, *, angle_deg: int, density: float, max_step: int,
               collision_model: str = "hybrid", dt_s: float | None = None,
               beam_diag_interval: int, background_diag_interval: int,
               beam_particles: int | None = None, random_seed: int = 42) -> Path:
    if dt_s is None:
        dt_s = PARAMETERS["dt_s"]
    if beam_particles is None:
        beam_particles = PARAMETERS["beam_particles"]
    case_dir = RUNS / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "inputs").write_text(
        input_text(
            angle_deg=angle_deg,
            density=density,
            max_step=max_step,
            dt_s=dt_s,
            collision_model=collision_model,
            beam_diag_interval=beam_diag_interval,
            background_diag_interval=background_diag_interval,
            beam_particles=beam_particles,
            random_seed=random_seed,
        ),
        encoding="utf-8",
    )
    script = case_dir / "run.sh"
    script.write_text(run_script(case_name), encoding="utf-8", newline="\n")
    script.chmod(0o755)
    manifest = {
        "case": case_name,
        "angle_deg": angle_deg,
        "background_density_m3": density,
        "background_ion_temperature_eV": PARAMETERS["background_ion_temperature_eV"],
        "background_electron_temperature_eV": PARAMETERS["background_electron_temperature_eV"],
        "max_step": max_step,
        "dt_s": dt_s,
        "end_time_s": max_step * dt_s,
        "beam_particles": beam_particles,
        "random_seed": random_seed,
        "collision_model": collision_model,
        "ion_coulomb_enabled": collision_model in {"hybrid", "ion_only"},
        "electron_stopping_enabled": collision_model in {"hybrid", "electron_only"},
        "model": "fixed-position kinetic proton reservoir with selectable collisions",
        "background_position_model": PARAMETERS["background_position_model"],
    }
    (case_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return case_dir


def main() -> None:
    RUNS.mkdir(exist_ok=True)
    smoke = PARAMETERS["smoke"]
    smoke_dir = write_case(
        "smoke_fixed_angle_035_n_1e17",
        angle_deg=smoke["angle_deg"],
        density=smoke["background_density_m3"],
        max_step=smoke["max_step"],
        beam_diag_interval=smoke["beam_diag_interval"],
        background_diag_interval=smoke["background_diag_interval"],
    )

    validation_spec = [
        ("validation_angle_035_n_1e15_hybrid", 1.0e15, "hybrid", PARAMETERS["dt_s"], 100),
        ("validation_angle_035_n_1e16_hybrid", 1.0e16, "hybrid", PARAMETERS["dt_s"], 100),
        ("validation_angle_035_n_1e17_ion_only", 1.0e17, "ion_only", PARAMETERS["dt_s"], 100),
        ("validation_angle_035_n_1e17_electron_only", 1.0e17, "electron_only", PARAMETERS["dt_s"], 100),
        ("validation_angle_035_n_1e17_half_dt", 1.0e17, "hybrid", 0.5 * PARAMETERS["dt_s"], 200),
    ]
    validation = [
        write_case(
            name,
            angle_deg=35,
            density=density,
            collision_model=model,
            dt_s=dt_s,
            max_step=max_step,
            beam_diag_interval=10,
            background_diag_interval=20,
        )
        for name, density, model, dt_s, max_step in validation_spec
    ]

    pilot = [
        write_case(
            f"pilot_angle_{angle:03d}_n_1e16_hybrid",
            angle_deg=angle,
            density=1.0e16,
            collision_model="hybrid",
            max_step=100,
            beam_diag_interval=10,
            background_diag_interval=20,
        )
        for angle in PARAMETERS["angles_deg"]
    ]
    controls = [
        write_case(
            f"control_angle_{angle:03d}_n_1e16_none",
            angle_deg=angle,
            density=1.0e16,
            collision_model="none",
            max_step=100,
            beam_diag_interval=10,
            background_diag_interval=20,
        )
        for angle in PARAMETERS["angles_deg"]
    ]
    horizon = write_case(
        "horizon_angle_035_n_1e16_hybrid",
        angle_deg=35,
        density=1.0e16,
        collision_model="hybrid",
        max_step=500,
        beam_diag_interval=10,
        background_diag_interval=100,
    )

    production = []
    for angle in PARAMETERS["angles_deg"]:
        for density in PARAMETERS["background_densities_m3"]:
            name = f"angle_{angle:03d}_n_{density_tag(density)}"
            production.append(
                write_case(
                    name,
                    angle_deg=angle,
                    density=density,
                    max_step=PARAMETERS["production_max_step"],
                    beam_diag_interval=PARAMETERS["beam_diag_interval"],
                    background_diag_interval=PARAMETERS["background_diag_interval"],
                )
            )

    smoke_launcher = RUNS / "run_smoke.sh"
    smoke_launcher.write_text(
        "#!/usr/bin/env bash\nset -e\ncd \"$(dirname \"$0\")/smoke_fixed_angle_035_n_1e17\"\n./run.sh\n",
        encoding="utf-8",
        newline="\n",
    )
    smoke_launcher.chmod(0o755)

    validation_launcher = RUNS / "run_validation.sh"
    validation_lines = ["#!/usr/bin/env bash", "set -e", 'root="$(cd "$(dirname "$0")" && pwd)"']
    for case in validation:
        validation_lines.append(f'(cd "$root/{case.name}" && ./run.sh)')
    validation_launcher.write_text(
        "\n".join(validation_lines) + "\n", encoding="utf-8", newline="\n"
    )
    validation_launcher.chmod(0o755)

    pilot_launcher = RUNS / "run_pilot.sh"
    pilot_lines = ["#!/usr/bin/env bash", "set -e", 'root="$(cd "$(dirname "$0")" && pwd)"']
    for case in pilot:
        pilot_lines.append(f'(cd "$root/{case.name}" && ./run.sh)')
    pilot_launcher.write_text(
        "\n".join(pilot_lines) + "\n", encoding="utf-8", newline="\n"
    )
    pilot_launcher.chmod(0o755)

    control_launcher = RUNS / "run_pilot_controls.sh"
    control_lines = ["#!/usr/bin/env bash", "set -e", 'root="$(cd "$(dirname "$0")" && pwd)"']
    for case in controls:
        control_lines.append(f'(cd "$root/{case.name}" && ./run.sh)')
    control_launcher.write_text(
        "\n".join(control_lines) + "\n", encoding="utf-8", newline="\n"
    )
    control_launcher.chmod(0o755)

    horizon_launcher = RUNS / "run_horizon_probe.sh"
    horizon_launcher.write_text(
        f'#!/usr/bin/env bash\nset -e\ncd "$(dirname "$0")/{horizon.name}"\n./run.sh\n',
        encoding="utf-8",
        newline="\n",
    )
    horizon_launcher.chmod(0o755)

    production_launcher = RUNS / "run_production.sh"
    production_lines = ["#!/usr/bin/env bash", "set -e", 'root="$(cd "$(dirname "$0")" && pwd)"']
    for case in production:
        production_lines.append(f'(cd "$root/{case.name}" && ./run.sh)')
    production_launcher.write_text(
        "\n".join(production_lines) + "\n", encoding="utf-8", newline="\n"
    )
    production_launcher.chmod(0o755)

    print(f"Generated smoke case: {smoke_dir.relative_to(ROOT)}")
    print(f"Generated {len(validation)} focused validation cases")
    print(f"Generated {len(pilot)} moderate-density angle-pilot cases")
    print(f"Generated {len(controls)} matched collision-off controls")
    print(f"Generated horizon probe: {horizon.relative_to(ROOT)}")
    print(f"Generated {len(production)} production cases below {RUNS.relative_to(ROOT)}")
    print("Production cases were generated but are not started automatically.")


if __name__ == "__main__":
    main()
