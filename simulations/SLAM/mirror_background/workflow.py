#!/usr/bin/env python3
"""Generate, run, and audit GPU background-plasma extensions of Lars's source.

Run with the existing WSL warpx-gpu Python. Runtime decks and large diagnostics
live in ignored runs/. No original input, field map, or WarpX source is modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
LARS = REPO / "simulations/Mirror/inputs_mirror.txt"
STL = LARS.parent / "cylinder_OGMirror.stl"
FIELD_DEFAULT = REPO.parent / "hackathon-beam-anisotropy/warpx/official_magnetic_mirror/example-femm-3d.h5"
WARPX_DEFAULT = Path("/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d")
QE, ME, MP, EPS0, C = 1.602176634e-19, 9.1093837139e-31, 1.67262192595e-27, 8.8541878188e-12, 299792458.0
MODELS = ("ion_background", "ion_coulomb", "hybrid", "kinetic", "kinetic_half_dt")


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def field_audit(path: Path) -> dict:
    """Check native openPMD HDF5 geometry; no assumed SI scale or array order."""
    with h5py.File(path, "r") as data:
        fields = []
        data.visititems(lambda name, obj: fields.append(obj.name) if isinstance(obj, h5py.Group)
                        and name.endswith("/B") and all(c in obj for c in ("x", "y", "z")) else None)
        if len(fields) != 1:
            raise ValueError(f"Expected one magnetic mesh, found {fields}")
        group = data[fields[0]]
        vectors = [np.asarray(group[c], dtype=float) * float(group[c].attrs.get("unitSI", 1.0)) for c in ("x", "y", "z")]
        if not all(np.isfinite(a).all() for a in vectors):
            raise ValueError("Non-finite magnetic field")
        if len({a.shape for a in vectors}) != 1:
            raise ValueError("Staggered component array shapes require a separate field audit")
        axes = [a.decode() if isinstance(a, bytes) else str(a) for a in group.attrs["axisLabels"]]
        spacing = np.asarray(group.attrs["gridSpacing"], float) * float(group.attrs.get("gridUnitSI", 1.0))
        offset = np.asarray(group.attrs["gridGlobalOffset"], float) * float(group.attrs.get("gridUnitSI", 1.0))
        position = np.asarray(group["x"].attrs.get("position", [0, 0, 0]), float)
        for component in ("y", "z"):
            if not np.array_equal(position, group[component].attrs.get("position", [0, 0, 0])):
                raise ValueError("Staggered component positions require a separate field audit")
        order = group.attrs.get("dataOrder", "C")
        if isinstance(order, bytes):
            order = order.decode()
        if order != "C" or sorted(axes) != ["x", "y", "z"]:
            raise ValueError("This workflow supports Cartesian C-order xyz-labelled fields only")
        low = offset + spacing * position
        high = low + spacing * (np.asarray(vectors[0].shape) - 1)
        bounds = {axis: [float(low[i]), float(high[i])] for i, axis in enumerate(axes)}
        for axis, required in zip(("x", "y", "z"), ((-1, 1), (-1, 1), (0, 5))):
            if bounds[axis][0] > required[0] + 1e-10 or bounds[axis][1] < required[1] - 1e-10:
                raise ValueError(f"Field does not cover simulation on {axis}: {bounds[axis]}")
        return {"path": str(path), "sha256": digest(path), "mesh": fields[0],
                "axis_labels": axes, "shape": list(vectors[0].shape), "bounds_m": bounds,
                "grid_spacing_m": spacing.tolist(),
                "Bmax_T": float(np.sqrt(sum(v*v for v in vectors)).max()),
                "identity": "Located legacy official mirror map; Lars's exact copy not independently confirmed",
                "periodicity": "Off-axis transverse B is not z-periodic; see BOUNDARY_STATUS.md"}


def source_parameters() -> dict[str, str]:
    """Retain the live continuous source parameters, excluding old commented beam."""
    values = {}
    for line in LARS.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    if values.get("protons.injection_style") != "NFluxPerCell":
        raise ValueError("Lars's source changed; inspect it before reusing this extension")
    return values


def collision(name: str, species: str) -> str:
    return f"""{name}.type = pairwisecoulomb
{name}.species = {species}
{name}.CoulombLog = 10
{name}.ndt_supercycle = 1
"""


def background(name: str, species: str, temperature: float, kinetic: bool, ppc: int) -> str:
    mass = ME if species == "electron" else MP
    ustd = math.sqrt(temperature * QE / mass) / C
    return f"""{name}.species_type = {species}
{name}.do_not_deposit = {0 if kinetic else 1}
{name}.do_not_push = {0 if kinetic else 1}
{name}.initialize_self_fields = 0
{name}.injection_style = NUniformPerCell
{name}.num_particles_per_cell_each_dim = {ppc} {ppc} {ppc}
{name}.profile = parse_density_function
{name}.density_function(x,y,z) = "n_bg*(x*x+y*y<0.75*0.75)*(z>0.25)*(z<4.75)"
{name}.momentum_distribution_type = gaussian
{name}.ux_th = {ustd:.17g}
{name}.uy_th = {ustd:.17g}
{name}.uz_th = {ustd:.17g}
"""


def generate(case: str, run: Path, field: Path, steps: int | None, ppc: int, wall: str) -> dict:
    if case not in MODELS or wall not in ("stl", "analytic-bore"):
        raise ValueError("Unknown case or wall model")
    kinetic = case.startswith("kinetic")
    audit = field_audit(field)
    dt_limit = min(0.1 * ME / (QE * audit["Bmax_T"]),
                   0.1 / math.sqrt(1e10 * QE**2 / (ME*EPS0)),
                   0.2 * 0.0625 / (4*math.sqrt(10*QE/ME)))
    dt = min(5e-11, dt_limit) if kinetic else min(1e-7, 0.1*MP/(QE*audit["Bmax_T"]))
    nsteps = steps if steps is not None else 200
    if case == "kinetic_half_dt":
        dt *= 0.5
        nsteps *= 2
    density = 1e10 if kinetic else 1e16
    debye_e = math.sqrt(EPS0*10/(density*QE))
    debye_combined = debye_e / math.sqrt(11)
    plasma_dt = math.sqrt(density*QE**2/(ME*EPS0))*dt
    gyro_dt = QE*audit["Bmax_T"]/ME*dt
    if kinetic and (debye_combined < 0.0625 or plasma_dt > 0.1 + 1e-12 or gyro_dt > 0.1 + 1e-12):
        raise ValueError("Kinetic resolution check failed")
    original = source_parameters()
    source_keys = [k for k in original if k.startswith("protons.") and k not in (
        "protons.do_not_deposit", "protons.initialize_self_fields")]
    source = "\n".join(f"{k} = {original[k]}" for k in source_keys)
    constants = "\n".join(f"{k} = {original[k]}" for k in (
        "my_constants.y0", "my_constants.z0", "my_constants.sigma_y", "my_constants.sigma_z"))
    names, operators = [], []
    if case != "ion_background":
        names.append("beam_ion")
        operators.append(collision("beam_ion", "protons background_ions"))
    if case == "hybrid":
        names.append("electron_drag")
        operators.append("""electron_drag.type = background_stopping
electron_drag.species = protons
electron_drag.background_type = electrons
electron_drag.background_density(x,y,z,t) = "n_bg*(x*x+y*y<0.75*0.75)*(z>0.25)*(z<4.75)"
electron_drag.background_mass = m_e
electron_drag.background_temperature = 10*q_e/kb
""")
    if kinetic:
        for name, pair in (("beam_electron", "protons electrons"), ("electron_ion", "electrons background_ions"),
                           ("ion_ion", "background_ions background_ions"), ("electron_electron", "electrons electrons"),
                           ("beam_beam", "protons protons")):
            names.append(name)
            operators.append(collision(name, pair))
    species = "protons background_ions" + (" electrons" if kinetic else "")
    diagnostic_interval = max(1, nsteps // 4)
    wall_input = (f'eb2.geom_type = stl\neb2.stl_file = "{STL}"\neb2.cover_multiple_cuts = 1'
                  if wall == "stl" else 'warpx.eb_implicit_function = "x*x+y*y-0.9*0.9"')
    text = f"""# Generated by workflow.py from Lars's continuous source. Case: {case}
# Explicit standard-periodic exploratory variant, NOT Periodic_ReflectParVel.
# Background starts away from the nonperiodic field seam. Short integration only.
# Reservoir cases omit deposition; kinetic cases enable self-consistent E.
geometry.dims = 3
geometry.prob_lo = -1 -1 0
geometry.prob_hi = 1 1 5
amr.n_cell = 32 32 80
amr.max_level = 0
amr.max_grid_size = 128
amr.blocking_factor = 8
max_step = {nsteps}
warpx.const_dt = {dt:.17g}
warpx.random_seed = 42
warpx.serialize_initial_conditions = 1
warpx.verbose = 1
warpx.use_filter = 0
warpx.do_electrostatic = labframe
warpx.grid_type = collocated
algo.particle_shape = 1
warpx.self_fields_required_precision = 1e-8
warpx.self_fields_absolute_tolerance = 1e-10
warpx.eb_potential(x,y,z,t) = "0.0"
{wall_input}
boundary.field_lo = pec pec periodic
boundary.field_hi = pec pec periodic
boundary.particle_lo = absorbing absorbing periodic
boundary.particle_hi = absorbing absorbing periodic
boundary.particle_eb = absorbing
my_constants.n_bg = {density:.17g}
{constants}
particles.species_names = {species}
protons.do_not_deposit = {0 if kinetic else 1}
protons.initialize_self_fields = 0
{source}

{background('background_ions', 'proton', 1.0, kinetic, ppc)}
"""
    if kinetic:
        text += "\n" + background("electrons", "electron", 10.0, True, ppc)
    if names:
        text += "\ncollisions.correct_energy_momentum = 0\ncollisions.collision_names = " + " ".join(names) + "\n" + "\n".join(operators)
    # Huge beam/background weight differences: record uncorrected operator choice.
    text += f"""
particles.B_ext_particle_init_style = read_from_file
particles.read_fields_from_path = "{field}"
diagnostics.diags_names = diag
diag.diag_type = Full
diag.format = openpmd
diag.openpmd_backend = h5
diag.file_prefix = diags/diag/openpmd
diag.intervals = {diagnostic_interval}
diag.fields_to_plot = Ex Ey Ez rho phi Bx By Bz eb_covered
diag.species = {species}
diag.write_species = 1
"""
    for name in species.split():
        text += f"diag.{name}.variables = x y z ux uy uz w\n"
    text += "\nwarpx.reduced_diags_names = number energy momentum field_energy field_max\nreduced_diags.path = diags/reduced/\nreduced_diags.precision = 18\n"
    for name, typ in (("number", "ParticleNumber"), ("energy", "ParticleEnergy"), ("momentum", "ParticleMomentum"),
                      ("field_energy", "FieldEnergy"), ("field_max", "FieldMaximum")):
        text += f"{name}.type = {typ}\n{name}.intervals = 1\n"
    (run / "inputs.txt").write_text(text, encoding="utf-8")
    metadata = {"case": case, "workflow_sha256": digest(Path(__file__)),
                "source_input": str(LARS), "source_sha256": digest(LARS),
                "stl": str(STL), "stl_sha256": digest(STL), "field": audit,
                "dt_s": dt, "max_step": nsteps, "horizon_s": dt*nsteps, "density_m3": density,
                "Ti_eV": 1, "Te_eV": 10, "ppc_each_dimension": ppc, "species": species.split(),
                "electron_mass_kg": ME if kinetic else None,
                "self_consistent_E": kinetic, "analytic_electron_drag": case == "hybrid",
                "collision_names": names, "pairwise_CoulombLog": 10,
                "collision_moment_correction": False,
                "background_spatial_profile": "r<0.75 m, 0.25<z<4.75 m; identical ion/electron loading",
                "boundary": "explicit standard periodic z; custom boundary unavailable",
                "wall": wall, "wall_potential_V": 0.0,
                "poisson_relative_tolerance": 1e-8, "poisson_absolute_tolerance_V_m2": 1e-10,
                "cell_size_m": 0.0625, "lambda_De_m": debye_e, "lambda_D_combined_m": debye_combined,
                "electron_omega_p_dt": plasma_dt, "electron_omega_c_max_dt": gyro_dt,
                "kinetic_resolution_checks_apply": kinetic,
                "limitations": ["Short integration smoke, not confinement or collision-rate validation",
                    "HDF5 seam is not periodic: do not use wrapped trajectories as physical evidence",
                    "Kinetic and reservoir densities differ; not a controlled physical comparison",
                    "Original tiny beam flux retained; source physical strength needs experimental specification",
                    "Unequal-weight pair collisions conserve statistically, not exactly per event",
                    "GPU random seed does not guarantee bitwise reproducibility"]}
    write_json(run / "metadata.json", metadata)
    return metadata


def read_reduced(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open() as stream:
        header = stream.readline()
    labels = re.findall(r"\[\d+\]([^\s]+)", header)
    return labels, np.loadtxt(path, ndmin=2)


def analyze(run: Path) -> dict:
    """Operational checks + honest diagnostics; not a claim of convergence."""
    metadata = json.loads((run / "metadata.json").read_text())
    execution = json.loads((run / "execution.json").read_text())
    result = {"case": metadata["case"], "run": str(run), "execution": execution,
              "horizon_s": metadata["horizon_s"], "checks": {}, "species": {}}
    checks = result["checks"]
    checks["exit_zero"] = execution["returncode"] == 0
    log = (run / "run.log").read_text(errors="replace")
    checks["cuda_initialized"] = "CUDA initialized with" in log
    if not checks["exit_zero"]:
        result["status"] = "failed"
        result["log_tail"] = log[-6000:]
        write_json(run / "summary.json", result)
        return result
    if metadata["collision_names"]:
        checks["binary_collision_dispatch_in_profile"] = bool(re.search(
            r"BinaryCollision::doCollisionsWithinTile::LoopOverCollisions\s+[1-9][0-9]*", log))
    if metadata["analytic_electron_drag"]:
        checks["analytic_stopping_dispatch_in_profile"] = bool(re.search(
            r"BackgroundStopping::doCollisions\(\)\s+[1-9][0-9]*", log))
    if metadata["self_consistent_E"]:
        checks["no_analytic_stopping_double_count"] = "BackgroundStopping::doCollisions()" not in log
    reduced = {}
    for name in ("number", "energy", "momentum", "field_energy", "field_max"):
        labels, values = read_reduced(run / f"diags/reduced/{name}.txt")
        reduced[name] = {label: values[:, i] for i, label in enumerate(labels)}
        checks[f"{name}_finite"] = bool(np.isfinite(values).all())
    checks["horizon_reached"] = bool(np.isclose(reduced["number"]["time(s)"][-1], metadata["horizon_s"], rtol=1e-8, atol=1e-18))
    from openpmd_viewer import OpenPMDTimeSeries
    series = OpenPMDTimeSeries(str(run / "diags/diag/openpmd"), backend="h5py")
    checks["first_snapshot_at_t0"] = bool(abs(float(series.t[0])) < 1e-18)
    result["output_snapshots"] = len(series.iterations)
    background_initial = {}
    for species in metadata["species"]:
        n = reduced["number"][f"{species}_macroparticles()"]
        energy = reduced["energy"][f"{species}_mean(J)"]
        weight = reduced["number"][f"{species}_weight()"]
        item = {"first_reduced_count": int(n[0]), "final_count": int(n[-1]),
                "final_physical_weight": float(weight[-1]), "final_mean_energy_eV": float(energy[-1]/QE),
                "initial_mean_energy_eV": float(energy[0]/QE)}
        checks[f"{species}_nonempty"] = bool(n[-1] > 0)
        # Full diagnostic t=0 backgrounds are useful for exact neutral loading.
        if species != "protons":
            first = series.get_particle(["x", "y", "z", "w"], species=species, iteration=series.iterations[0])
            background_initial[species] = first
        clearance = float("inf")
        finite = True
        for iteration in series.iterations:
            arrays = series.get_particle(["x", "y", "z", "ux", "uy", "uz", "w"], species=species, iteration=iteration)
            finite = finite and all(np.isfinite(a).all() for a in arrays)
            if arrays[2].size:
                clearance = min(clearance, float(np.min(arrays[2])), float(5-np.max(arrays[2])))
        checks[f"{species}_snapshots_finite"] = bool(finite)
        checks[f"{species}_away_from_z_seam_at_snapshots"] = bool(clearance > 0.05)
        item["min_snapshot_z_seam_clearance_m"] = clearance if math.isfinite(clearance) else None
        result["species"][species] = item
        if species != "protons":
            checks[f"{species}_no_losses_short_run"] = bool(n[-1] == n[0])
            target = 15 if species == "electrons" else 1.5
            checks[f"{species}_initial_thermal_energy"] = bool(abs(energy[0]/QE/target-1) < 0.03)
    if metadata["self_consistent_E"]:
        ion = background_initial["background_ions"]
        electron = background_initial["electrons"]
        # Sort to ignore diagnostic particle storage order across GPU tiles.
        def ordered(arrays):
            xyz = np.column_stack(arrays)
            return xyz[np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0]))]
        a, b = ordered(ion), ordered(electron)
        checks["initial_background_exactly_neutral_sampling"] = bool(a.shape == b.shape and np.allclose(a, b, rtol=1e-13, atol=1e-15))
        result["initial_background_net_charge_C"] = float(QE*(np.sum(ion[3])-np.sum(electron[3])))
        initial_ex, _ = series.get_field(field="E", coord="x", iteration=series.iterations[0])
        final_ex, _ = series.get_field(field="E", coord="x", iteration=series.iterations[-1])
        result["initial_max_abs_Ex_V_m"] = float(np.max(np.abs(initial_ex)))
        result["final_max_abs_Ex_V_m"] = float(np.max(np.abs(final_ex)))
        # Tolerance here is V/m, independent of the Poisson residual's V/m^2.
        checks["electric_response_above_1e-8_V_m"] = bool(np.max(np.abs(final_ex)) > max(1e-8, 100*np.max(np.abs(initial_ex))))
        result["resolution"] = {key: metadata[key] for key in ("lambda_De_m", "lambda_D_combined_m", "cell_size_m", "electron_omega_p_dt", "electron_omega_c_max_dt")}
        kinetic_energy = reduced["energy"]["total(J)"]
        electric_keys = [key for key in reduced["field_energy"] if key.startswith("E_lev")]
        if len(electric_keys) == 1:
            total = kinetic_energy + reduced["field_energy"][electric_keys[0]]
            result["remaining_particle_plus_E_energy_relative_change"] = float(total[-1]/total[0]-1)
            result["energy_accounting_note"] = "Remaining energy only; injected/escaped energy not subtracted. This short run has tiny beam weight."
    result["status"] = "integration_smoke_pass" if all(checks.values()) else "needs_review"
    result["not_validated"] = metadata["limitations"]
    write_json(run / "summary.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "analyze"))
    parser.add_argument("--case", choices=MODELS, default="kinetic")
    parser.add_argument("--field", type=Path, default=FIELD_DEFAULT)
    parser.add_argument("--warpx", type=Path, default=WARPX_DEFAULT)
    parser.add_argument("--boundary", choices=("standard-periodic",), help="Required acknowledgment: custom boundary is not implemented")
    parser.add_argument("--steps", type=int)
    parser.add_argument("--ppc", type=int, default=2)
    parser.add_argument("--wall", choices=("stl", "analytic-bore"), default="analytic-bore")
    parser.add_argument("--tag", default="")
    parser.add_argument("--run-dir", type=Path)
    args = parser.parse_args()
    if args.action == "analyze":
        if not args.run_dir:
            parser.error("--run-dir is required for analyze")
        print(json.dumps(analyze(args.run_dir.resolve()), indent=2))
        return
    if args.boundary != "standard-periodic":
        parser.error("Select --boundary standard-periodic explicitly; read BOUNDARY_STATUS.md")
    if args.steps is not None and args.steps <= 0 or args.ppc <= 0:
        parser.error("Steps and ppc must be positive")
    if not re.fullmatch(r"[A-Za-z0-9_-]*", args.tag):
        parser.error("Tag must contain only letters, digits, underscores, or hyphens")
    stamp = time.strftime("%Y%m%d-%H%M%S")
    run = HERE / "runs" / (args.case + "_" + (args.tag or stamp))
    run.mkdir(parents=True, exist_ok=False)
    metadata = generate(args.case, run, args.field.resolve(strict=True), args.steps, args.ppc, args.wall)
    print(f"Prepared {run}; {metadata['max_step']} steps, dt={metadata['dt_s']:.4g} s", flush=True)
    if args.action == "prepare":
        return
    executable = args.warpx.resolve(strict=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    command = [str(executable), "inputs.txt"]
    start = time.monotonic()
    with (run / "run.log").open("w") as log:
        process = subprocess.run(command, cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT)
    execution = {"returncode": process.returncode, "wall_seconds": time.monotonic()-start,
                 "command": command, "binary_sha256": digest(executable),
                 "source_commit": subprocess.run(["git", "-C", str(executable.parents[2]), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()}
    write_json(run / "execution.json", execution)
    result = analyze(run)
    print(json.dumps(result, indent=2), flush=True)
    if result["status"] != "integration_smoke_pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
