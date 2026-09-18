"""Separate true hybrid-PIC validation case; never modifies legacy cases.

The analytic periodic mirror is a solver/diagnostic fixture, NOT a SLAM field
or a reproduction of Lars's FEMM mirror. Source parameters are reused verbatim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

from diagnostics import Grid, regions, write_json

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parents[1] / "Mirror" / "inputs_mirror.txt"
BINARY = Path("/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d")
QE, MP, C = 1.602176634e-19, 1.67262192595e-27, 299792458.0


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def build_input(steps=20, dt=1e-9, ppc=3, closure="isothermal", collisions=False,
                diagnostic_every=None, adaptive_fields=False):
    max_dt = 5e-9 if adaptive_fields else 1e-9
    if steps < 1 or not math.isfinite(dt) or dt <= 0 or dt > max_dt or ppc < 2:
        raise ValueError(
            "Validation case: positive steps, ppc >= 2 per axis, "
            f"0 < dt <= {max_dt*1e9:g} ns; dt above 1 ns requires adaptive fields")
    if diagnostic_every is not None and (
            isinstance(diagnostic_every, bool) or
            not isinstance(diagnostic_every, int) or diagnostic_every <= 0):
        raise ValueError("diagnostic_every must be a positive integer number of steps")
    if closure not in ("isothermal", "polytropic", "energy_transport"):
        raise ValueError("Choose isothermal, polytropic, or energy_transport")
    energy_equation = closure == "energy_transport"
    cadence = max(1, steps // 4) if diagnostic_every is None else diagnostic_every
    values = {}
    for line in SOURCE.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            values[k.strip()] = v.strip()
    if values.get("protons.injection_style") != "NFluxPerCell":
        raise ValueError("Lars's live source changed; review before generating")
    if (values.get("protons.flux_normal_axis") != "x" or
            float(values.get("protons.surface_flux_pos", "nan")) != -.2 or
            "my_constants.y0" not in values):
        raise ValueError("This fixture expects the OLD mirror source, not the new SLAM source; review migration")
    source = "\n".join(f"{k} = {v}" for k, v in values.items()
                       if k.startswith("protons.") and k not in
                       ("protons.do_not_deposit", "protons.initialize_self_fields"))
    constants = "\n".join(f"{k} = {values[k]}" for k in
                          ("my_constants.y0", "my_constants.z0",
                           "my_constants.sigma_y", "my_constants.sigma_z"))
    gamma = 1.0 if closure == "isothermal" else 5 / 3
    closure_comment = (
        "# Electron temperature evolves by entropy transport and compression/expansion.\n"
        "# No collisional or Joule heating; no heat conduction or physical electron wall loss.\n"
        "# Temperature changes alone must NOT be attributed to beam heating.\n"
        "# See ELECTRON_MODEL_AUDIT.md before enabling any electron-ion relaxation."
        if energy_equation else
        "# Isothermal/polytropic closure, NOT an electron-heating calculation.")
    adaptive_controls = (
        "# RKF45 adapts FIELD substeps only; particle/fluid-advection dt stays const_dt.\n"
        "hybrid_pic_model.use_rkf45 = 1\n"
        "hybrid_pic_model.substep_rtol = 1.e-5\n"
        "hybrid_pic_model.substep_atol = 1.e-12\n"
        if adaptive_fields else "")
    fluid_velocity_diagnostics = (
        "# Electron FLUID velocity in m/s, not a kinetic-electron velocity distribution.\n"
        "# Explicit registry components bypass this build's generic vector-lookup bug.\n"
        'diag.additional_fields_to_plot = "hybrid_electron_velocity_fp[dir=x]" '
        '"hybrid_electron_velocity_fp[dir=y]" "hybrid_electron_velocity_fp[dir=z]"\n'
        if energy_equation else "")
    u_th = math.sqrt(QE / MP) / C
    text = f'''# TRUE HYBRID-PIC VALIDATION FIXTURE -- NOT SLAM / NOT THE FEMM MAP
# Electrons: inertialess quasi-neutral fluid. Protons: moving kinetic particles.
{closure_comment}
geometry.dims = 3
geometry.prob_lo = -1 -1 0
geometry.prob_hi = 1 1 5
amr.n_cell = 32 32 80
amr.max_level = 0
amr.max_grid_size = 128
amr.blocking_factor = 8
max_step = {steps}
warpx.const_dt = {dt:.17g}
warpx.random_seed = 42
warpx.serialize_initial_conditions = 1
warpx.verbose = 1
warpx.grid_type = collocated
warpx.use_filter = 1
algo.particle_shape = 1
algo.current_deposition = direct
algo.maxwell_solver = hybrid
# Field subcycling remains necessary for Hall/whistler stability.
hybrid_pic_model.substeps = 64
{adaptive_controls}# Reference Te in eV; evolved when the electron-energy equation is enabled.
hybrid_pic_model.elec_temp = 10
hybrid_pic_model.gamma = {gamma:.17g}
hybrid_pic_model.n0_ref = 1.e16
hybrid_pic_model.n_floor = 1.e15
hybrid_pic_model.holmstrom_vacuum_region = 1
hybrid_pic_model.plasma_resistivity(rho,J,t) = "0.0"
hybrid_pic_model.plasma_hyper_resistivity(rho,B) = "0.0"
hybrid_pic_model.solve_electron_energy_equation = {int(energy_equation)}
# Resistivity is zero and Joule heating is explicitly off.
hybrid_pic_model.include_joule_heating = 0
# Deliberately omit the electron-ion relaxation-rate parameter: installed
# drag/thermal exchange does not close the independent beam drift-energy budget.

# Divergence-free applied mirror: B = curl A, smoothly periodic in z.
# B_axis minimum 0.001487675 T at z=2.5; maximum/minimum = 3.10951921.
# The external A split makes the field available to BOTH fluid and ions.
my_constants.b_min = 0.001487675
my_constants.mirror_ratio = 3.10951921
my_constants.b_avg = b_min*(mirror_ratio+1)/2
my_constants.a_mirror = (mirror_ratio-1)/(mirror_ratio+1)
my_constants.kz = 2*pi/5
hybrid_pic_model.add_external_fields = 1
external_vector_potential.fields = mirror
external_vector_potential.do_diva_cleaning = 0
external_vector_potential.mirror.Ax_external_grid_function(x,y,z) = "-0.5*y*b_avg*(1+a_mirror*cos(kz*z))"
external_vector_potential.mirror.Ay_external_grid_function(x,y,z) = "0.5*x*b_avg*(1+a_mirror*cos(kz*z))"
external_vector_potential.mirror.Az_external_grid_function(x,y,z) = "0.0"
external_vector_potential.mirror.A_time_external_function(t) = "1.0"

# Analytic validation bore, not Lars's SLAM vessel. No Debye sheath model.
warpx.eb_implicit_function = "x*x+y*y-0.9*0.9"
boundary.field_lo = pec pec periodic
boundary.field_hi = pec pec periodic
boundary.particle_lo = absorbing absorbing periodic
boundary.particle_hi = absorbing absorbing periodic
boundary.particle_eb = absorbing
my_constants.n_bg = 1.e16
{constants}
particles.species_names = protons background_ions
protons.do_not_deposit = 0
protons.initialize_self_fields = 0
{source}

# Uniform 1-eV protons inside r<0.75 m; not a force-balanced equilibrium.
# Full periodic z extent removes artificial axial ends of the background.
background_ions.species_type = proton
background_ions.do_not_deposit = 0
background_ions.do_not_push = 0
background_ions.initialize_self_fields = 0
background_ions.injection_style = NUniformPerCell
background_ions.num_particles_per_cell_each_dim = {ppc} {ppc} {ppc}
background_ions.profile = parse_density_function
background_ions.density_function(x,y,z) = "n_bg*(x*x+y*y<0.75*0.75)"
background_ions.momentum_distribution_type = gaussian
background_ions.ux_th = {u_th:.17g}
background_ions.uy_th = {u_th:.17g}
background_ions.uz_th = {u_th:.17g}

# Keep both populations separate; physical weights are used in analysis.
diagnostics.diags_names = diag
diag.diag_type = Full
diag.format = openpmd
diag.openpmd_backend = h5
diag.file_prefix = diags/diag/openpmd
diag.intervals = {cadence}
# Te is K (convert with kB/qe); Pe is Pa. Separate species and physical weights.
diag.fields_to_plot = Ex Ey Ez Bx By Bz rho Te Pe eb_covered
{fluid_velocity_diagnostics}diag.species = protons background_ions
diag.write_species = 1
diag.protons.variables = x y z ux uy uz w
diag.background_ions.variables = x y z ux uy uz w
'''
    if collisions:
        text += '''
# Optional ion-ion Coulomb operators; no electron drag/energy exchange.
# Extreme beam/background weight ratio: exact energy correction disabled.
collisions.correct_energy_momentum = 0
collisions.collision_names = beam_ion ion_ion
beam_ion.type = pairwisecoulomb
beam_ion.species = protons background_ions
beam_ion.CoulombLog = 10
ion_ion.type = pairwisecoulomb
ion_ion.species = background_ions background_ions
ion_ion.CoulombLog = 10
'''
    manifest = dict(model="true kinetic-ion / inertialess-fluid-electron hybrid PIC",
                    geometry="analytic periodic mirror validation fixture, NOT SLAM or FEMM",
                    closure=closure, gamma=gamma, Te0_eV=10, Ti0_eV=1,
                    n0_m3=1e16, n_floor_m3=1e15, dt_s=dt, steps=steps,
                    horizon_s=dt*steps, field_substeps=64, ppc_per_axis=ppc,
                    ion_collisions=collisions, electron_heating=False,
                    electron_heating_definition="No electron-ion collisional or Joule heat sources; compression can change Te",
                    electron_energy_equation=energy_equation,
                    electron_temperature_evolves=closure != "isothermal",
                    electron_velocity_output_fields=(
                        [f"hybrid_electron_velocity_fp[dir={axis}]" for axis in "xyz"]
                        if energy_equation else []),
                    adaptive_fields=bool(adaptive_fields),
                    field_substep_rtol=1e-5 if adaptive_fields else None,
                    field_substep_atol=1e-12 if adaptive_fields else None,
                    diagnostic_every_steps=cadence,
                    source_path=str(SOURCE), source_sha256=sha(SOURCE),
                    grid=Grid().as_dict(), regions=regions(Grid()),
                    beam_flux_m2_s=1e5, beam_area_m2=0.04,
                    beam_rate_physical_s=4000,
                    limits=["No electron velocity PDF or kinetic electron physics",
                            "No electron-ion drag or collisional energy transfer enabled",
                            "No equilibrium initialization, no heating attribution",
                            "Periodic recirculation is not end-loss confinement",
                            ("Electron entropy transport and compression; no conduction or physical electron wall heat loss"
                             if energy_equation else
                             "Algebraic pressure closure, not a transported electron energy equation"),
                            "No beam-heating attribution from a changing temperature alone",
                            "Adaptive field steps do not relax particle or electron-fluid advection timestep limits",
                            "No AMR; density floor and vacuum switch require sensitivity tests"])
    return text, manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", type=Path, help="New output directory; existing paths are refused")
    p.add_argument("--steps", type=int, default=20)
    p.add_argument("--dt", type=float, default=1e-9)
    p.add_argument("--ppc", type=int, default=3, help="markers per axis (3 -> 27 per cell)")
    p.add_argument("--closure", choices=("isothermal", "polytropic", "energy_transport"), default="isothermal")
    p.add_argument("--diagnostic-every", type=int, default=None,
                   help="Positive output cadence in steps; default max(1, steps // 4)")
    p.add_argument("--adaptive-fields", action="store_true",
                   help="RKF45 field substeps (rtol=1e-5, atol=1e-12); permits dt up to 5 ns")
    p.add_argument("--ion-collisions", action="store_true")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--binary", type=Path, default=BINARY)
    args = p.parse_args()
    deck, manifest = build_input(args.steps, args.dt, args.ppc, args.closure,
                                 args.ion_collisions, args.diagnostic_every,
                                 args.adaptive_fields)
    args.run.mkdir(parents=True, exist_ok=False)
    (args.run / "inputs.txt").write_text(deck, encoding="utf-8")
    if args.execute:
        manifest["binary"] = str(args.binary.resolve())
        manifest["binary_sha256"] = sha(args.binary)
    write_json(args.run / "manifest.json", manifest)
    print(f"Prepared {args.run.resolve()}", flush=True)
    if args.execute:
        start = time.monotonic()
        with (args.run / "warpx.log").open("w") as log:
            result = subprocess.run([str(args.binary), "inputs.txt"], cwd=args.run,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    env={**os.environ, "OMP_NUM_THREADS": "4"})
        write_json(args.run / "execution.json", dict(returncode=result.returncode,
                   wall_seconds=time.monotonic()-start))
        print(f"WarpX exit: {result.returncode}", flush=True)
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
