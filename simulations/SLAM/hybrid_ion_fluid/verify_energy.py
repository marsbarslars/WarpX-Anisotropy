"""Small closed-box verification of the installed hybrid electron-energy operator.

This unmagnetized, periodic uniform box is a NUMERICAL TEST, not the physical
mirror and not a lab collision-rate prediction. Its constant Qei rate is
deliberately artificial. A thermal test and a counterstream energy-accounting
test distinguish temperature relaxation from bulk-drift energy transfer.
No simulation runs unless --execute is explicitly supplied to prepare.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np

QE, KB, MP, C = 1.602176634e-19, 1.380649e-23, 1.67262192595e-27, 299792458.
MU0 = 1.25663706212e-6
BINARY = Path("/home/ryanv/src/warpx/build_gpu_py/bin/warpx.3d")


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def build_input(case="thermal", steps=100, dt=1e-9, ppc=4, cells=8, nu=1e5):
    if case not in ("thermal", "counterstream"):
        raise ValueError("case must be thermal or counterstream")
    if steps < 10 or cells not in (8, 16) or not 3 <= ppc <= 8:
        raise ValueError("Need >=10 steps, 8 or 16 cells/axis and 3..8 markers/axis")
    if not math.isfinite(dt) or not 0 < dt <= 1e-9:
        raise ValueError("Verification timestep must be positive and <=1 ns")
    if not math.isfinite(nu) or nu < 0 or nu*dt > .001:
        raise ValueError("Rate must be nonnegative and nu*dt <=0.001")
    names = ["ions"] if case == "thermal" else ["plus", "minus"]
    te, ti, n0, side, gamma = (10. if case == "thermal" else 1.), 1., 1e16, .5, 5/3
    v0 = 0. if case == "thermal" else 1e5
    spread = math.sqrt(QE*ti/MP)/C
    text = f'''# UNIFORM PERIODIC NUMERICAL VERIFICATION -- NOT A MIRROR OR LAB PREDICTION
# Artificial constant nu={nu:g}/s; no fitted physical collision rate.
geometry.dims = 3
geometry.prob_lo = 0 0 0
geometry.prob_hi = {side} {side} {side}
amr.n_cell = {cells} {cells} {cells}
amr.max_level = 0
amr.max_grid_size = 32
amr.blocking_factor = 8
max_step = {steps}
warpx.const_dt = {dt:.17g}
warpx.random_seed = 173
warpx.serialize_initial_conditions = 1
warpx.verbose = 1
warpx.grid_type = collocated
warpx.use_filter = 1
algo.particle_shape = 1
algo.current_deposition = direct
algo.maxwell_solver = hybrid
hybrid_pic_model.substeps = 10
hybrid_pic_model.elec_temp = {te}
hybrid_pic_model.gamma = {gamma:.17g}
hybrid_pic_model.n0_ref = {n0:.17g}
hybrid_pic_model.n_floor = {n0*.05:.17g}
hybrid_pic_model.holmstrom_vacuum_region = 0
hybrid_pic_model.solve_electron_energy_equation = 1
hybrid_pic_model.include_joule_heating = 0
hybrid_pic_model.plasma_resistivity(rho,J,t) = "0.0"
hybrid_pic_model.plasma_hyper_resistivity(rho,B) = "0.0"
hybrid_pic_model.electron_ion_relaxation_rate(rho,Te,Ti,t) = "{nu:.17g}"
# No imposed B, source, wall, beam injection, loss or ion-ion operator.
boundary.field_lo = periodic periodic periodic
boundary.field_hi = periodic periodic periodic
boundary.particle_lo = periodic periodic periodic
boundary.particle_hi = periodic periodic periodic
particles.species_names = {" ".join(names)}
'''
    for index, name in enumerate(names):
        bulk = v0*(1 if index == 0 else -1)/C
        text += f'''
{name}.species_type = proton
{name}.do_not_deposit = 0
{name}.do_not_push = 0
{name}.do_temperature_deposition = 1
{name}.initialize_self_fields = 0
{name}.injection_style = NUniformPerCell
{name}.num_particles_per_cell_each_dim = {ppc} {ppc} {ppc}
{name}.profile = constant
{name}.density = {n0/len(names):.17g}
{name}.momentum_distribution_type = gaussian
{name}.ux_m = {bulk:.17g}
{name}.uy_m = 0
{name}.uz_m = 0
{name}.ux_th = {spread:.17g}
{name}.uy_th = {spread:.17g}
{name}.uz_th = {spread:.17g}
'''
    text += f'''
diagnostics.diags_names = diag
diag.diag_type = Full
diag.format = openpmd
diag.openpmd_backend = h5
diag.file_prefix = diags/diag/openpmd
# Include the first initialized snapshot; iteration 0 Pe is not initialized.
diag.intervals = 1:1,{max(1, steps//10)},{steps}:{steps}
diag.fields_to_plot = Ex Ey Ez Bx By Bz rho Te Pe {" ".join("T_"+name for name in names)}
diag.species = {" ".join(names)}
diag.write_species = 1
'''
    for name in names:
        text += f"diag.{name}.variables = x y z ux uy uz w\n"
    manifest = dict(case=case, purpose="NUMERICAL VERIFICATION ONLY, not the physical mirror",
                    artificial_rate_s_inv=nu, physical_collision_rate=False,
                    steps=steps, dt_s=dt, horizon_s=steps*dt, cells_per_axis=cells,
                    side_m=side, ppc_per_axis=ppc, species=names, n0_m3=n0,
                    gamma=gamma, Te0_eV=te, Ti0_eV=ti, drift_speed_m_s=v0,
                    reference="Uniform one-proton-species d(Te-Ti)/dt=-4*nu*(Te-Ti)",
                    limitations=["Finite-marker stochastic and deposition errors remain",
                                 "Zero magnetic field: not a confinement test",
                                 "Counterstream runs audit drift-energy accounting; thermal formula is not their complete solution",
                                 "No claim of exact discrete energy conservation"])
    return text, manifest


def thermal_reference(times, te0, ti0, nu):
    times = np.asarray(times, float)
    if np.any(~np.isfinite(times)) or np.any(times < 0) or nu < 0:
        raise ValueError("Need finite nonnegative elapsed times and rate")
    mean = .5*(te0+ti0)
    difference = (te0-ti0)*np.exp(-4*nu*times)
    return mean+.5*difference, mean-.5*difference


def particle_budget(u, weights):
    u, weights = np.asarray(u, float), np.asarray(weights, float)
    if u.shape != (len(weights), 3) or not len(weights):
        raise ValueError("Need nonempty Nx3 normalized momenta and N weights")
    if not np.isfinite(u).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Need finite momenta and positive finite weights")
    u2 = np.einsum("ij,ij->i", u, u)
    lorentz = np.sqrt(1+u2)
    velocity = C*u/lorentz[:, None]
    if np.linalg.norm(velocity, axis=1).max() > .01*C:
        raise ValueError("Temperature moments are nonrelativistic")
    number = float(weights.sum())
    mean = np.sum(weights[:, None]*velocity, axis=0)/number
    thermal = .5*MP*float(weights @ np.sum((velocity-mean)**2, axis=1))
    return dict(physical_count=number, marker_count=len(weights),
                K_J=MP*C*C*float(weights @ (u2/(lorentz+1))),
                bulk_K_J=.5*MP*number*float(mean@mean),
                Ti_eV=thermal/(1.5*number*QE), mean_vx_m_s=float(mean[0]))


def relative_change(values):
    values = np.asarray(values, float)
    if len(values) == 0 or not np.isfinite(values).all() or values[0] <= 0:
        raise ValueError("Need finite energy history with positive reference")
    return (values-values[0])/values[0]


def analyze(run):
    from openpmd_viewer import OpenPMDTimeSeries
    run = Path(run)
    manifest = json.loads((run/"manifest.json").read_text())
    execution = json.loads((run/"execution.json").read_text())
    if execution["returncode"]:
        raise ValueError("Refusing analysis of an unsuccessful run")
    series = OpenPMDTimeSeries(str(run/"diags/diag/openpmd"), backend="h5py")
    if int(series.iterations[-1]) != manifest["steps"]:
        raise ValueError("Missing final requested iteration")
    if set(series.avail_species) != set(manifest["species"]):
        raise ValueError("Unexpected particle populations")
    volume, rows, skipped = manifest["side_m"]**3, [], []
    for iteration, t in zip(series.iterations, series.t):
        def field(name, coord=None):
            result, _ = series.get_field(field=name, coord=coord, iteration=int(iteration))
            array = np.asarray(result, float)
            if array.shape != (manifest["cells_per_axis"],)*3 or not np.isfinite(array).all():
                raise ValueError(f"Unexpected or nonfinite field {name}")
            return array
        pe, te, rho = field("Pe"), field("Te")*KB/QE, field("rho")
        if iteration == 0:
            skipped.append(dict(step=0, reason="Pre-hybrid initialization; not a fluid-energy reference"))
            continue
        if np.any(pe <= 0) or np.any(te <= 0) or np.any(rho <= 0):
            raise ValueError("Uniform verification plasma must have positive Pe,Te,rho")
        row = dict(step=int(iteration), time_s=float(t),
                   Te_eV=float(np.sum(rho*te)/np.sum(rho)),
                   Pe_mean_Pa=float(pe.mean()), rho_mean_C_m3=float(rho.mean()),
                   density_relative_std=float(rho.std()/rho.mean()),
                   electron_U_J=float(pe.mean()*volume/(manifest["gamma"]-1)))
        populations = []
        for species in manifest["species"]:
            ux, uy, uz, weights = series.get_particle(["ux", "uy", "uz", "w"],
                                                    species=species, iteration=int(iteration))
            budget = particle_budget(np.column_stack((ux, uy, uz)), weights)
            populations.append(budget)
            row.update({f"{species}_{key}": value for key, value in budget.items()})
            # Native grid Ti is a different estimator (deposition + local drift).
            row[f"{species}_deposited_Ti_mean_eV"] = float(field("T_"+species).mean())
        row["ion_K_J"] = sum(p["K_J"] for p in populations)
        row["ion_bulk_K_J"] = sum(p["bulk_K_J"] for p in populations)
        row["ion_Ti_eV"] = sum(p["Ti_eV"]*p["physical_count"] for p in populations)/sum(
            p["physical_count"] for p in populations)
        row["ion_plus_electron_J"] = row["ion_K_J"]+row["electron_U_J"]
        row["magnetic_U_J"] = volume/(2*MU0)*sum(float(np.mean(field("B", c)**2)) for c in "xyz")
        row["E_max_V_m"] = max(float(np.max(np.abs(field("E", c)))) for c in "xyz")
        row["ion_electron_magnetic_J"] = row["ion_plus_electron_J"]+row["magnetic_U_J"]
        rows.append(row)
    if len(rows) < 3:
        raise ValueError("Need at least three initialized snapshots")
    ts = np.array([r["time_s"] for r in rows])
    te = np.array([r["Te_eV"] for r in rows])
    ti = np.array([r["ion_Ti_eV"] for r in rows])
    energy = np.array([r["ion_plus_electron_J"] for r in rows])
    residual = relative_change(energy)
    for row, value in zip(rows, residual):
        row["energy_relative_change"] = float(value)
    summary = dict(case=manifest["case"], artificial_rate_s_inv=manifest["artificial_rate_s_inv"],
                   reference_step=rows[0]["step"], reference_time_s=rows[0]["time_s"],
                   skipped=skipped, final_time_s=float(ts[-1]),
                   signed_final_energy_relative_change=float(residual[-1]),
                   max_abs_energy_relative_change=float(np.max(np.abs(residual))),
                   bulk_energy_change_J=float(rows[-1]["ion_bulk_K_J"]-rows[0]["ion_bulk_K_J"]),
                   electron_energy_change_J=float(rows[-1]["electron_U_J"]-rows[0]["electron_U_J"]),
                   magnetic_energy_max_J=max(r["magnetic_U_J"] for r in rows),
                   marker_counts_constant=all(all(r[f"{s}_marker_count"] == rows[0][f"{s}_marker_count"]
                       for s in manifest["species"]) for r in rows),
                   no_exact_conservation_claim=True,
                   interpretation="Signed closed-box residual; finite-marker noise and operator splitting remain. Compare zero-rate control before assigning cause.")
    reference = None
    if manifest["case"] == "thermal":
        reference = thermal_reference(ts-ts[0], te[0], ti[0], manifest["artificial_rate_s_inv"])
        delta0 = te[0]-ti[0]
        misfit = ((te-ti)-(reference[0]-reference[1]))/abs(delta0)
        summary["max_normalized_temperature_difference_misfit"] = float(np.max(np.abs(misfit)))
        for row, trefe, trefi, error in zip(rows, *reference, misfit):
            row.update(Te_reference_eV=float(trefe), Ti_reference_eV=float(trefi),
                       temperature_difference_misfit=float(error))
    with (run/"verification_history.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    write_json(run/"verification_summary.json", summary)
    plot(run, manifest, rows, reference)
    print(json.dumps(summary, indent=2))
    return summary


def plot(run, manifest, rows, reference):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.7))
    ns = np.array([r["time_s"] for r in rows])*1e9
    ax[0].plot(ns, [r["Te_eV"] for r in rows], "o-", ms=3, label="Fluid electrons")
    ax[0].plot(ns, [r["ion_Ti_eV"] for r in rows], "s-", ms=3, label="Ions: drift-subtracted")
    if reference is not None:
        ax[0].plot(ns, reference[0], "--", color="black", label="Uniform thermal reference")
        ax[0].plot(ns, reference[1], "--", color="black")
        ax[1].plot(ns, [100*r["temperature_difference_misfit"] for r in rows], "o-")
        ax[1].set(ylabel="Difference misfit / initial difference (%)", title="Thermal-reference misfit")
    else:
        base = rows[0]["ion_plus_electron_J"]
        ax[1].plot(ns, [(r["ion_bulk_K_J"]-rows[0]["ion_bulk_K_J"])/base*100 for r in rows], label="Ion drift energy")
        ax[1].plot(ns, [(r["electron_U_J"]-rows[0]["electron_U_J"])/base*100 for r in rows], label="Electron internal energy")
        ax[1].set(ylabel="Energy change / initial total (%)", title="Where does drift energy go?")
        ax[1].legend(fontsize=8)
    ax[0].set(ylabel="Temperature (eV)", title="Measured temperatures")
    ax[0].legend(fontsize=8)
    ax[2].plot(ns, [r["energy_relative_change"]*100 for r in rows], "o-", color="C3")
    ax[2].set(ylabel="(Ion K + electron U − initial) / initial (%)", title="Signed closed-box energy residual")
    for axes in ax:
        axes.set_xlabel("Time (ns)"); axes.grid(alpha=.22)
    ax[1].axhline(0, color="gray", lw=.8); ax[2].axhline(0, color="gray", lw=.8)
    fig.suptitle(f"{manifest['case'].title()} verification | artificial ν = {manifest['artificial_rate_s_inv']:g} s⁻¹ | NOT a mirror prediction")
    fig.text(.5, .005, "Reference starts at first initialized snapshot. Finite-marker noise and splitting remain; no exact conservation claim.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .94))
    fig.savefig(run/"energy_verification.png", dpi=180)
    fig.savefig(run/"energy_verification.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("run", type=Path)
    prepare.add_argument("--case", choices=("thermal", "counterstream"), default="thermal")
    prepare.add_argument("--steps", type=int, default=100)
    prepare.add_argument("--dt", type=float, default=1e-9)
    prepare.add_argument("--ppc", type=int, default=4)
    prepare.add_argument("--cells", type=int, default=8)
    prepare.add_argument("--nu", type=float, default=1e5, help="ARTIFICIAL verification rate; use 0 for paired control")
    prepare.add_argument("--execute", action="store_true")
    prepare.add_argument("--binary", type=Path, default=BINARY)
    analyze_command = commands.add_parser("analyze")
    analyze_command.add_argument("run", type=Path)
    args = parser.parse_args()
    if args.action == "analyze":
        analyze(args.run)
        return
    deck, manifest = build_input(args.case, args.steps, args.dt, args.ppc, args.cells, args.nu)
    args.run.mkdir(parents=True, exist_ok=False)
    (args.run/"inputs.txt").write_text(deck, encoding="utf-8")
    if args.execute:
        with args.binary.open("rb") as binary:
            manifest["binary_sha256"] = hashlib.file_digest(binary, "sha256").hexdigest()
        manifest["binary"] = str(args.binary.resolve())
    write_json(args.run/"manifest.json", manifest)
    print(f"Prepared {args.run.resolve()}", flush=True)
    if args.execute:
        start = time.monotonic()
        with (args.run/"warpx.log").open("w") as log:
            result = subprocess.run([str(args.binary), "inputs.txt"], cwd=args.run,
                                    stdout=log, stderr=subprocess.STDOUT,
                                    env={**os.environ, "OMP_NUM_THREADS": "4"})
        write_json(args.run/"execution.json", dict(returncode=result.returncode,
                   wall_seconds=time.monotonic()-start))
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
