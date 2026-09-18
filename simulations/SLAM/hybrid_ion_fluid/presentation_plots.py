"""Slide-ready figures from measured hybrid diagnostics; never interpolate data.

Usage: python presentation_plots.py RUN [--reference MATCHED_RUN] [--output DIR]
Creates PNG/SVG figures and numerical comparison records, not a slide deck.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

QE = 1.602176634e-19
REGIONS = ("center", "near_injection")
REGION_LABELS = {"center": "Mirror center", "near_injection": "Near injection"}
MOMENTS = ("T_parallel_eV", "T_perp_eV", "T_scalar_eV", "P_parallel_Pa",
           "P_perp_Pa", "P_scalar_Pa", "density_m3")
NONPHYSICS = {"dt_s", "steps", "horizon_s", "source_path", "binary", "limits",
              "diagnostic_interval", "diagnostic_interval_steps", "output_interval_steps",
              "diagnostic_every", "diagnostic_every_steps", "diagnostics_every", "diagnostics_every_steps",
              "output_every", "output_interval", "write_interval"}
REQUIRED_PHYSICS = ("model", "geometry", "closure", "source_sha256", "grid",
                    "n0_m3", "n_floor_m3", "Te0_eV", "Ti0_eV", "ppc_per_axis",
                    "field_substeps", "ion_collisions")


def close_time(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-18)


def read_rows(path):
    rows = []
    with Path(path).open(newline="", encoding="utf-8") as stream:
        for original in csv.DictReader(stream):
            row = {}
            for key, value in original.items():
                if value in ("", None):
                    row[key] = None
                elif key in ("region", "species"):
                    row[key] = value
                elif value in ("True", "False"):
                    row[key] = value == "True"
                else:
                    row[key] = float(value)
                    if not math.isfinite(row[key]):
                        raise ValueError(f"Nonfinite {key} in {path}")
            rows.append(row)
    if not rows:
        raise ValueError(f"Empty diagnostic table: {path}")
    return rows


def selection(rows, region, species=None, initialized=False):
    found = [r for r in rows if r["region"] == region and
             (species is None or r.get("species") == species) and
             (not initialized or r.get("initialized_pressure", r["step"] > 0))]
    found.sort(key=lambda r: r["time_s"])
    if any(close_time(a["time_s"], b["time_s"]) for a, b in zip(found, found[1:])):
        raise ValueError(f"Duplicate diagnostic times for {region}/{species}")
    return found


def load_run(path):
    path = Path(path).resolve()
    data = {"path": path,
            "manifest": json.loads((path/"manifest.json").read_text()),
            "audit": json.loads((path/"diagnostic_audit.json").read_text()),
            "ions": read_rows(path/"ion_history.csv"),
            "fluid": read_rows(path/"electron_fluid_history.csv")}
    if data["audit"].get("finite") is not True:
        raise ValueError("A successful finite-data diagnostic audit is required")
    if (path/"execution.json").exists():
        if json.loads((path/"execution.json").read_text())["returncode"] != 0:
            raise ValueError("Refusing presentation figures from a failed run")
    for region in REGIONS:
        for table, species in (("ions", "background_ions"), ("fluid", None)):
            rows = selection(data[table], region, species)
            if not rows or not close_time(rows[-1]["time_s"], data["manifest"]["horizon_s"]):
                raise ValueError(f"Missing requested endpoint: {table}/{region}")
    return data


def compare_background(run, reference):
    """Return signed scaled differences at coincident samples only.

    All physics/provenance fields must agree; numerical timestep/step count and
    output cadence may differ. The reference must be identified, not called truth.
    """
    a, b = run["manifest"], reference["manifest"]
    missing = [key for key in REQUIRED_PHYSICS if key not in a or key not in b]
    if missing:
        raise ValueError(f"Missing comparison physics metadata: {missing}")
    differences = [key for key in sorted((set(a) | set(b))-NONPHYSICS) if a.get(key) != b.get(key)]
    if differences:
        raise ValueError(f"Different scenarios cannot be a timestep misfit: {differences}")
    if not close_time(a["horizon_s"], b["horizon_s"]):
        raise ValueError("Reference must have the same physical horizon")
    if not (a["dt_s"] > 0 and b["dt_s"] > 0):
        raise ValueError("Both comparison timesteps must be positive")
    records, unmatched = [], {}
    for region in REGIONS:
        rr = selection(run["ions"], region, "background_ions")
        ref = selection(reference["ions"], region, "background_ions")
        if not rr or not ref or not close_time(rr[0]["time_s"], 0) or not close_time(ref[0]["time_s"], 0):
            raise ValueError("Both runs need initial background moments")
        if not close_time(rr[-1]["time_s"], a["horizon_s"]) or not close_time(ref[-1]["time_s"], b["horizon_s"]):
            raise ValueError("Comparison data do not reach matching endpoints")
        matched = []
        for row in rr:
            candidates = [r for r in ref if close_time(row["time_s"], r["time_s"])]
            if not candidates:
                continue
            r = candidates[0]
            matched.append(row["time_s"])
            for key in MOMENTS:
                if row.get(key) is None or r.get(key) is None:
                    raise ValueError(f"Undefined background moment {key} at a comparison time")
                if key.startswith("T_"):
                    scale, definition = b["Ti0_eV"], "reference input Ti0 (eV)"
                elif key.startswith("P_"):
                    scale, definition = ref[0][key], f"reference measured initial {key}"
                else:
                    scale, definition = b["n0_m3"], "reference input n0 (1/m^3)"
                if scale is None or not math.isfinite(scale) or scale <= 0:
                    raise ValueError(f"Positive initial normalization required for {key}")
                records.append({"region": region, "time_s": row["time_s"], "quantity": key,
                                "run_value": row[key], "reference_value": r[key],
                                "scale": scale, "scale_definition": definition,
                                "signed_scaled_difference": (row[key]-r[key])/scale})
        if len(matched) < 2 or not close_time(matched[-1], a["horizon_s"]):
            raise ValueError("At least initial and final coincident samples are required")
        unmatched[region] = {"matched_times_s": matched,
                             "unmatched_run_samples": len(rr)-len(matched),
                             "unmatched_reference_samples": len(ref)-len(matched)}
    fluid_records, fluid_sampling = compare_fluid(run, reference)
    return {"definition": "(run - reference) / documented initial reference scale",
            "interpolation": False, "reference_is_truth": False,
            "run_dt_s": a["dt_s"], "reference_dt_s": b["dt_s"],
            "horizon_s": a["horizon_s"], "sampling": unmatched, "records": records,
            "max_abs_scaled_difference": max(abs(r["signed_scaled_difference"]) for r in records),
            "fluid_records": fluid_records, "fluid_sampling": fluid_sampling,
            "max_abs_fluid_scaled_difference": max(abs(r["signed_scaled_difference"]) for r in fluid_records),
            "interpretation": "Numerical comparison of short-run background-ion and electron-fluid moments; not experimental accuracy, equilibrium, or heating validation"}


def compare_fluid(run, reference):
    """Compare initialized, coincident fluid samples after scenario validation.

    Called by compare_background, which checks model identity and equal horizons.
    The reference is a numerical comparator, not a known exact solution.
    """
    manifest = reference["manifest"]
    records, sampling = [], {}
    for region in REGIONS:
        rows = selection(run["fluid"], region, initialized=True)
        ref = selection(reference["fluid"], region, initialized=True)
        if not rows or not ref or not close_time(rows[-1]["time_s"], manifest["horizon_s"]) or not close_time(ref[-1]["time_s"], manifest["horizon_s"]):
            raise ValueError("Both fluid comparisons require initialized final endpoint data")
        keys = ["Te_eV", "Pe_Pa", "ne_grid_m3"]
        if all(r.get("Te_effective_eV") is not None for r in rows+ref):
            keys.insert(1, "Te_effective_eV")
        matched = []
        for row in rows:
            candidates = [r for r in ref if close_time(row["time_s"], r["time_s"])]
            if not candidates:
                continue
            other = candidates[0]
            matched.append(row["time_s"])
            for key in keys:
                if row.get(key) is None or other.get(key) is None:
                    raise ValueError(f"Undefined initialized electron-fluid quantity: {key}")
                if key.startswith("Te_"):
                    scale, definition = manifest["Te0_eV"], "reference input Te0 (eV)"
                elif key == "Pe_Pa":
                    scale = manifest["n0_m3"]*QE*manifest["Te0_eV"]
                    definition = "reference n0 * elementary charge * Te0 (Pa)"
                else:
                    scale, definition = manifest["n0_m3"], "reference input n0 (1/m^3)"
                if not math.isfinite(scale) or scale <= 0:
                    raise ValueError(f"Positive fluid comparison normalization required: {key}")
                records.append({"region": region, "time_s": row["time_s"], "quantity": key,
                                "run_value": row[key], "reference_value": other[key],
                                "scale": scale, "scale_definition": definition,
                                "signed_scaled_difference": (row[key]-other[key])/scale})
        if not matched or not close_time(matched[-1], manifest["horizon_s"]):
            raise ValueError("Coincident initialized fluid endpoint samples are required")
        sampling[region] = {"matched_times_s": matched,
                            "unmatched_run_samples": len(rows)-len(matched),
                            "unmatched_reference_samples": len(ref)-len(matched),
                            "initial_uninitialized_rows_omitted": True,
                            "effective_temperature_compared": "Te_effective_eV" in keys}
    return records, sampling


def time_axis(horizon):
    return (1e9, "ns") if horizon < 1e-6 else (1e6, "microseconds")


def load_pdf_payload(run):
    """Read saved radial velocity PDFs; missing legacy artifacts may be skipped.

    Values are already normalized by each entire population's physical weight,
    not by the in-range fraction. Do not renormalize or fit the distributions.
    """
    path = run["path"]/"velocity_distributions.npz"
    if not path.exists():
        return None
    panels = []
    with np.load(path, allow_pickle=False) as archive:
        par = np.asarray(archive["parallel_edges_m_s"], dtype=float)
        perp = np.asarray(archive["perp_edges_m_s"], dtype=float)
        if any(e.ndim != 1 or len(e) < 2 or not np.isfinite(e).all() or np.any(np.diff(e) <= 0)
               for e in (par, perp)) or perp[0] < 0:
            raise ValueError("Invalid velocity histogram edges")
        for region, species in (("center", "background_ions"), ("near_injection", "protons")):
            rows = selection(run["ions"], region, species)
            if not rows:
                raise ValueError(f"Missing ion history for {region}/{species} PDF")
            last = rows[-1]
            key = f"{int(last['step'])}_{region}_{species}"
            if key not in archive:
                raise ValueError(f"Missing final velocity PDF: {key}")
            H = np.asarray(archive[key], dtype=float)
            if H.shape != (len(par)-1, len(perp)-1) or not np.isfinite(H).all() or np.any(H < 0):
                raise ValueError(f"Malformed saved velocity PDF: {key}")
            captured = float(np.sum(H*np.diff(par)[:, None]*np.diff(perp)[None, :]))
            if captured > 1+1e-9:
                raise ValueError("Saved probability distribution exceeds unit integral")
            expected = last.get("histogram_weight_fraction")
            if expected is not None and not math.isclose(captured, expected, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError("Saved PDF integral disagrees with recorded captured-weight fraction")
            if last.get("macro_count", 0) == 0 and H.any():
                raise ValueError("Nonzero PDF for an empty particle population")
            panels.append({"region": region, "species": species, "H": H.copy(),
                           "history": last, "captured_weight_fraction": captured})
    central_beam = selection(run["ions"], "center", "protons")
    return {"parallel_edges_m_s": par, "perp_edges_m_s": perp, "panels": panels,
            "central_beam_macro_count": central_beam[-1]["macro_count"] if central_beam else None}


def save_figure(fig, output, name, title, footer):
    fig.suptitle(title, x=.06, ha="left", fontsize=22, weight="bold")
    fig.text(.06, .025, footer, fontsize=11, color="#465063")
    fig.subplots_adjust(left=.07, right=.91, bottom=.18, top=.80, wspace=.48)
    for suffix in ("png", "svg"):
        fig.savefig(output/f"{name}.{suffix}", dpi=180, facecolor="white")


def make_plots(run, output, comparison=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                         "axes.titlesize": 16, "axes.labelsize": 13,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.formatter.useoffset": False, "lines.linewidth": 2.3,
                         "svg.fonttype": "none", "legend.fontsize": 10})
    factor, unit = time_axis(run["manifest"]["horizon_s"])
    m = run["manifest"]
    subtitle = f"{m['closure']} electrons | {m['horizon_s']*factor:g} {unit} | analytic periodic mirror test"
    files = []

    def base():
        return plt.subplots(1, 2, figsize=(14.22, 8))

    def decorate(ax, region):
        ax.set(title=REGION_LABELS[region], xlabel=f"Time ({unit})")
        ax.grid(alpha=.18)

    fig, axes = base()
    for ax, region in zip(axes, REGIONS):
        ions = selection(run["ions"], region, "background_ions")
        fluid = selection(run["fluid"], region, initialized=True)
        ts = [r["time_s"]*factor for r in ions]
        for key, label, color in (("T_parallel_eV", r"Ion $T_\parallel$", "#1864ab"),
                                   ("T_perp_eV", r"Ion $T_\perp$", "#2b8a3e")):
            ax.plot(ts, [r[key] for r in ions], label=label, color=color)
        ax.set_ylabel("Background ion temperature (eV)")
        ax.set_ylim(0, max(m["Ti0_eV"]*1.15, max(r[k] for r in ions for k in ("T_parallel_eV", "T_perp_eV"))*1.12))
        twin = ax.twinx()
        twin.spines["right"].set_visible(True)
        ft = [r["time_s"]*factor for r in fluid]
        twin.plot(ft, [r["Te_eV"] for r in fluid], color="#d9480f", label=r"Fluid $T_e$: volume mean")
        effective = [r.get("Te_effective_eV") for r in fluid]
        if effective and all(v is not None for v in effective):
            twin.plot(ft, effective, color="#a61e4d", linestyle="--", label=r"Fluid $T_e$: $\langle P_e\rangle/(e\langle n_e\rangle)$")
        twin.set_ylabel("Electron-fluid temperature (eV)", color="#d9480f")
        te_values = [r["Te_eV"] for r in fluid]
        if effective and all(v is not None for v in effective):
            te_values += effective
        twin.set_ylim(0, max(m["Te0_eV"]*1.15, max(te_values, default=m["Te0_eV"])*1.12))
        handles, labels = ax.get_legend_handles_labels()
        h2, l2 = twin.get_legend_handles_labels()
        ax.legend(handles+h2, labels+l2, loc="center left", frameon=False)
        decorate(ax, region)
    save_figure(fig, output, "01_temperature", "Ion and electron-fluid temperatures", subtitle+"\nSeparate vertical scales; ion flow removed cell by cell. Startup changes do not establish beam heating.")
    plt.close(fig); files.append("01_temperature")

    fig, axes = base()
    for ax, region in zip(axes, REGIONS):
        ions = selection(run["ions"], region, "background_ions")
        fluid = selection(run["fluid"], region, initialized=True)
        for key, label, color in (("P_parallel_Pa", r"Ion $P_\parallel$", "#1864ab"),
                                   ("P_perp_Pa", r"Ion $P_\perp$", "#2b8a3e")):
            ax.plot([r["time_s"]*factor for r in ions], [r[key]*1e3 for r in ions], label=label, color=color)
        ax.plot([r["time_s"]*factor for r in fluid], [r["Pe_Pa"]*1e3 for r in fluid], color="#d9480f", label=r"Fluid $P_e$")
        ax.set_ylabel("Pressure (mPa)"); ax.set_ylim(bottom=0)
        twin = ax.twinx(); twin.spines["right"].set_visible(True)
        twin.plot([r["time_s"]*factor for r in ions], [r["density_m3"]/m["n0_m3"] for r in ions], "--", color="#343a40", label=r"Ion $n/n_0$: cell counts")
        twin.plot([r["time_s"]*factor for r in fluid], [r["ne_grid_m3"]/m["n0_m3"] for r in fluid], ":", color="#868e96", label=r"Electron $n/n_0$: grid")
        twin.set_ylabel(r"Density / input $n_0$")
        density_values = [r["density_m3"]/m["n0_m3"] for r in ions]+[r["ne_grid_m3"]/m["n0_m3"] for r in fluid]
        twin.set_ylim(0, max(1.1, max(density_values)*1.12))
        handles, labels = ax.get_legend_handles_labels(); h2, l2 = twin.get_legend_handles_labels()
        ax.legend(handles+h2, labels+l2, loc="center left", frameon=False)
        decorate(ax, region)
    save_figure(fig, output, "02_pressure_density", "Background pressure and density", subtitle+"\nFluid quantities omit uninitialized t=0 pressure. Grid deposition and nearest-cell counts are different estimators.")
    plt.close(fig); files.append("02_pressure_density")

    fig, axes = base()
    for species, color, label in (("background_ions", "#1864ab", "Background ions"), ("protons", "#d9480f", "Injected protons")):
        snapshots = run["audit"]["snapshots"]
        times = [s["time_s"]*factor for s in snapshots]
        for ax, key in zip(axes, ("markers", "physical_count")):
            ax.plot(times, [s["populations"][species][key] for s in snapshots], label=label, color=color, marker="o", markersize=4)
    for ax, label, key in zip(axes, ("Computational markers", "Represented physical particles"), ("markers", "physical_count")):
        positives = [s["populations"][species][key] for s in run["audit"]["snapshots"]
                     for species in ("background_ions", "protons") if s["populations"][species][key] > 0]
        threshold = min(positives)/10 if positives else 1.
        ax.set(xlabel=f"Time ({unit})", ylabel=label, title=label)
        ax.set_yscale("symlog", linthresh=threshold)
        ax.set_ylim(bottom=0)
        ax.yaxis.get_major_locator().set_params(numticks=9)
        if key == "physical_count":
            ax.set_yticks([value for value in ax.get_yticks() if value == 0 or value >= threshold*100])
        ax.grid(alpha=.18); ax.legend(frameon=False)
    final_beam = run["audit"]["snapshots"][-1]["populations"]["protons"]
    axes[1].text(.04, .23, f"Final beam: {final_beam['physical_count']:.4g} physical protons\nrepresented by {final_beam['markers']:,} markers",
                 transform=axes[1].transAxes, fontsize=11, bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9})
    save_figure(fig, output, "03_particle_inventory", "Particle counts are not interchangeable", subtitle+"\nWhole-domain inventories; symlog axes include zero. Periodic recirculation is not particle loss; no steady-state claim.")
    plt.close(fig); files.append("03_particle_inventory")

    payload = load_pdf_payload(run)
    if payload is not None:
        fig, axes = base()
        par = payload["parallel_edges_m_s"]/1e3
        perp = payload["perp_edges_m_s"]/1e3
        for ax, panel in zip(axes, payload["panels"]):
            H = panel["H"]*1e6  # (m/s)^-2 -> (km/s)^-2, preserving probability.
            row = panel["history"]
            population = "Background ions" if panel["species"] == "background_ions" else "Injected protons"
            ax.set(title=f"{population}: {REGION_LABELS[panel['region']].lower()}",
                   xlabel=r"$v_\parallel$ (km/s)", ylabel=r"$|v_\perp|$ (km/s)",
                   xlim=(par[0], par[-1]), ylim=(perp[0], perp[-1]))
            if H.any():
                image = ax.pcolormesh(par, perp, H.T, shading="flat", vmin=0, cmap="magma")
                colorbar = fig.colorbar(image, ax=ax, pad=.025, fraction=.045)
                colorbar.set_label(r"Probability density / (km/s)$^2$", fontsize=11)
            else:
                ax.set_facecolor("#f1f3f5")
                text = "No particles in this region" if not row["macro_count"] else "No probability inside displayed velocity range"
                ax.text(.5, .5, text, ha="center", wrap=True, transform=ax.transAxes)
            low_weight = row.get("low_neff_weight_fraction")
            extra = f"\nWeight in low-count cells: {100*low_weight:.1f}%" if low_weight is not None else ""
            ax.text(.03, .96, f"{int(row['macro_count']):,} markers | captured weight: {100*panel['captured_weight_fraction']:.1f}%"+extra,
                    va="top", transform=ax.transAxes, fontsize=10,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9})
        central_note = "No beam particles at the center." if payload["central_beam_macro_count"] == 0 else "Beam and background remain separate populations."
        save_figure(fig, output, "05_velocity_space", "Measured ion velocity distributions", subtitle+
                    "\nLocal total B; lab-frame velocities; separate color scales. Radial PDFs, not Maxwellian fits."
                    "\n"+central_note+" Injected ions form a directed/nonthermal population; sparse cells are flagged.")
        # Leave space for the right colorbar's tick values and vertical label.
        fig.subplots_adjust(left=.07, right=.88, wspace=.50)
        for suffix in ("png", "svg"):
            fig.savefig(output/f"05_velocity_space.{suffix}", dpi=180, facecolor="white")
        plt.close(fig); files.append("05_velocity_space")

    if comparison:
        fig, axes = base()
        for ax, region in zip(axes, REGIONS):
            for key, label, style in (("T_parallel_eV", r"$T_\parallel/T_{i0}$", "-"),
                                      ("T_perp_eV", r"$T_\perp/T_{i0}$", "-"),
                                      ("P_scalar_Pa", r"$P_i/P_{i,ref}(0)$", "--"),
                                      ("density_m3", r"$n_i/n_0$", ":")):
                records = [r for r in comparison["records"] if r["region"] == region and r["quantity"] == key]
                ax.plot([r["time_s"]*factor for r in records], [r["signed_scaled_difference"] for r in records], linestyle=style, label=label, marker="o", markersize=4)
            ax.axhline(0, color="#868e96", linewidth=.8)
            ax.set_ylabel("(Run - reference) / initial scale")
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=False)
            ax.legend(frameon=False, loc="best"); decorate(ax, region)
        d1, d2 = comparison["run_dt_s"]*1e9, comparison["reference_dt_s"]*1e9
        save_figure(fig, output, "04_numerical_misfit", "Timestep comparison: background-ion moments", f"Run dt={d1:g} ns; reference dt={d2:g} ns | Coincident saved samples only; no interpolation.\nThis is numerical agreement, not experimental accuracy, steady state, or a validated heating prediction.")
        plt.close(fig); files.append("04_numerical_misfit")
        fig, axes = base()
        for ax, region in zip(axes, REGIONS):
            for key, label, style in (("Te_eV", r"Volume-mean $T_e/T_{e0}$", "-"),
                                      ("Te_effective_eV", r"Effective $T_e/T_{e0}$", "--"),
                                      ("Pe_Pa", r"$P_e/(n_0 e T_{e0})$", "-"),
                                      ("ne_grid_m3", r"$n_e/n_0$", ":")):
                records = [r for r in comparison["fluid_records"] if r["region"] == region and r["quantity"] == key]
                if not records:
                    continue
                ax.plot([r["time_s"]*factor for r in records], [r["signed_scaled_difference"] for r in records],
                        linestyle=style, label=label, marker="o", markersize=4)
            ax.axhline(0, color="#868e96", linewidth=.8)
            ax.set_ylabel("(Run - reference) / initial scale")
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useOffset=False)
            ax.legend(frameon=False, loc="best"); decorate(ax, region)
        save_figure(fig, output, "06_fluid_numerical_misfit", "Timestep comparison: electron-fluid response",
                    f"Run dt={d1:g} ns; reference dt={d2:g} ns | Initialized, coincident saved samples only; no interpolation.\nFluid timestep agreement is a numerical check, not proof of beam heating or a steady state.")
        plt.close(fig); files.append("06_fluid_numerical_misfit")
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run = load_run(args.run)
    reference = load_run(args.reference) if args.reference else None
    comparison = compare_background(run, reference) if reference else None
    output = args.output or args.run/"presentation"
    output.mkdir(parents=True, exist_ok=True)
    files = make_plots(run, output, comparison)
    if comparison:
        (output/"numerical_misfit.json").write_text(json.dumps(comparison, indent=2, allow_nan=False)+"\n")
        with (output/"numerical_misfit.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(comparison["records"][0]))
            writer.writeheader(); writer.writerows(comparison["records"])
        with (output/"fluid_misfit.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(comparison["fluid_records"][0]))
            writer.writeheader(); writer.writerows(comparison["fluid_records"])
    provenance = {"run": str(run["path"]), "reference": str(reference["path"]) if reference else None,
                  "figures": [f"{name}.{suffix}" for name in files for suffix in ("png", "svg")],
                  "geometry": run["manifest"]["geometry"], "closure": run["manifest"]["closure"],
                  "horizon_s": run["manifest"]["horizon_s"], "interpolated_or_smoothed": False,
                  "fluid_t0_pressure_excluded": True, "electron_velocity_pdf": "Not defined by a fluid model"}
    provenance["velocity_pdf_status"] = "Plotted measured saved ion PDFs" if "05_velocity_space" in files else "Legacy run has no velocity_distributions.npz; skipped"
    (output/"figure_manifest.json").write_text(json.dumps(provenance, indent=2)+"\n")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
