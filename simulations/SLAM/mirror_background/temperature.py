#!/usr/bin/env python3
"""Measure finite-region temperature, excluding mean flow, from openPMD particles.

This is a nonrelativistic second-moment temperature, not a Maxwellian fit or a
claim of thermal equilibrium. WarpX/openpmd-viewer ux is p/(mc), not velocity.
No particle is resampled, no curve is smoothed, and no missing step is invented.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from workflow import C, ME, MP, QE, digest, write_json


def velocity_from_u(u: np.ndarray) -> np.ndarray:
    """Convert normalized momentum gamma*v/c to physical velocity in m/s."""
    return C * u / np.sqrt(1 + np.sum(u * u, axis=1))[:, None]


def thermal_moments(v: np.ndarray, w: np.ndarray, mass: float) -> dict:
    """Population covariance weighted by physical particle multiplicity.

    Use centered velocities to avoid subtracting two large nearly equal raw
    moments. The denominator is sum(w), not the number of computational markers.
    Empty regions have undefined temperature (None/blank), never zero temperature.
    """
    v, w = np.asarray(v, float), np.asarray(w, float)
    if v.shape != (len(w), 3) or mass <= 0:
        raise ValueError("Expected N-by-3 velocities, N weights, and positive mass")
    if not np.isfinite(v).all() or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError("Nonfinite particles or nonpositive physical weights")
    result = {"macro_count": len(w), "physical_weight": float(w.sum()),
              "effective_count": 0.0, "T_eV": None, "mean_kinetic_energy_eV": None,
              "bulk_energy_eV": None}
    result.update({f"{key}{axis}{suffix}": None for axis in "xyz"
                   for key, suffix in (("T", "_eV"), ("v", "_mean_m_s"))})
    if not len(w):
        return result
    if np.max(np.linalg.norm(v, axis=1)) > 0.05 * C:
        raise ValueError("Nonrelativistic temperature diagnostic: v exceeds 0.05c")
    mean = np.average(v, axis=0, weights=w)
    variance = np.average((v - mean) ** 2, axis=0, weights=w)
    components = mass * variance / QE
    bulk = mass * float(mean @ mean) / (2 * QE)
    result.update({"effective_count": float(w.sum() ** 2 / (w @ w)),
                   "T_eV": float(components.mean()), "bulk_energy_eV": bulk,
                   "mean_kinetic_energy_eV": float(1.5 * components.mean() + bulk)})
    for i, axis in enumerate("xyz"):
        result[f"T{axis}_eV"] = float(components[i])
        result[f"v{axis}_mean_m_s"] = float(mean[i])
    return result


def region_mask(x: np.ndarray, y: np.ndarray, z: np.ndarray, probe: dict) -> np.ndarray:
    return ((x*x + y*y < probe["radius_m"] ** 2)
            & (z > probe["z_min_m"]) & (z < probe["z_max_m"]))


def snapshot_moments(series, iteration: int, species: str, probe: dict) -> tuple[dict, bool]:
    arrays = series.get_particle(["x", "y", "z", "ux", "uy", "uz", "w"],
                                 species=species, iteration=iteration)
    if not all(np.isfinite(a).all() for a in arrays):
        raise ValueError("Nonfinite particle data")
    keep = region_mask(*arrays[:3], probe)
    v = velocity_from_u(np.column_stack(arrays[3:6])[keep])
    moments = thermal_moments(v, arrays[6][keep], ME if species == "electrons" else MP)
    return moments, bool(keep.all())


def plot_temperature(rows: list[dict], metadata: dict, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    backgrounds = [s for s in metadata["species"] if s != "protons"]
    fig, axes = plt.subplots(len(backgrounds), 2, figsize=(11.5, 3.1*len(backgrounds)+1.4),
                             squeeze=False, gridspec_kw={"width_ratios": [1.5, 1]})
    for i, species in enumerate(backgrounds):
        selected = [r for r in rows if r["species"] == species]
        t = np.array([r["time_s"] for r in selected]) * 1e9
        temperatures = np.array([r["T_eV"] if r["T_eV"] is not None else np.nan for r in selected])
        counts = np.array([r["macro_count"] for r in selected])
        target = metadata["Te_eV"] if species == "electrons" else metadata["Ti_eV"]
        label = "Electrons" if species == "electrons" else "Background ions"
        color = "#D17427" if species == "electrons" else "#237FA4"
        ax = axes[i, 0]
        ax.plot(t, temperatures, color=color, lw=2, label="Measured central temperature")
        ax.axhline(target, color="#65707B", ls="--", lw=1.2, label="Input temperature")
        ax.set(title=label, ylabel="Temperature [eV]", xlabel="Time [ns]")
        # Keep at least a +/-10% scale: do not exaggerate startup fluctuations.
        finite = temperatures[np.isfinite(temperatures)]
        low = min(0.90*target, float(finite.min())*0.99) if finite.size else 0.90*target
        high = max(1.10*target, float(finite.max())*1.01) if finite.size else 1.10*target
        ax.set_ylim(low, high)
        ax.legend(fontsize=8, loc="upper right")
        count_ax = axes[i, 1]
        count_ax.plot(t, counts, color=color, lw=1.8)
        count_ax.set(title="Computational particles in sampling region",
                     ylabel="Macroparticle count", xlabel="Time [ns]")
        count_ax.set_ylim(0, max(1, float(counts.max())*1.2))
        for a in (ax, count_ax):
            a.grid(alpha=0.18)
            a.spines[["top", "right"]].set_visible(False)
            formatter = ScalarFormatter(useOffset=False)
            formatter.set_scientific(False)
            a.yaxis.set_major_formatter(formatter)
    probe = metadata["temperature_probe"]
    fig.suptitle("Mirror-center temperature history", fontsize=18, x=0.07, ha="left", y=0.97)
    fig.text(0.07, 0.92, f"r < {probe['radius_m']} m | {probe['z_min_m']} < z < {probe['z_max_m']} m | "
             f"physical-weighted velocity spread, mean flow removed", fontsize=10, color="#44515E")
    fig.text(0.07, 0.035, "Startup / diagnostic test, not evidence of beam heating. "
             "Particle exchange and finite sampling can change local temperature.\n"
             "No smoothing or interpolated samples. Dashed lines are initialization targets, not measured heating thresholds.",
             fontsize=9, color="#44515E")
    fig.subplots_adjust(left=0.08, right=0.97, top=0.83, bottom=0.17, hspace=0.5, wspace=0.30)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def analyze_temperature(run: Path) -> dict:
    from openpmd_viewer import OpenPMDTimeSeries

    metadata = json.loads((run / "metadata.json").read_text())
    probe = metadata.get("temperature_probe")
    if not probe:
        raise ValueError("This run has no central probe. Re-run with --temperature-every 1.")
    series = OpenPMDTimeSeries(str(run / probe["path"]), backend="h5py")
    full = OpenPMDTimeSeries(str(run / "diags/diag/openpmd"), backend="h5py")
    rows, checks, by_key = [], {}, {}
    expected = set(range(0, metadata["max_step"]+1, probe["interval_steps"]))
    expected.add(metadata["max_step"])  # FullDiagnostics writes final snapshot.
    checks["expected_iterations"] = set(map(int, series.iterations)) == expected
    checks["timestamps_match_steps"] = bool(np.allclose(series.t,
        np.asarray(series.iterations)*metadata["dt_s"], rtol=1e-10, atol=1e-18))
    all_inside = True
    for iteration, time_s in zip(series.iterations, series.t):
        for species in probe["species"]:
            moments, inside = snapshot_moments(series, int(iteration), species, probe)
            all_inside &= inside
            row = {"step": int(iteration), "time_s": float(time_s), "species": species, **moments}
            rows.append(row)
            by_key[(int(iteration), species)] = moments
    checks["probe_particles_inside_region"] = all_inside
    # Independently select the same region from the ordinary unfiltered snapshots.
    compared = 0
    for iteration in full.iterations:
        for species in probe["species"]:
            key = (int(iteration), species)
            if key not in by_key:
                continue
            reference, _ = snapshot_moments(full, int(iteration), species, probe)
            actual = by_key[key]
            match = all((actual[k] is None and reference[k] is None) if actual[k] is None or reference[k] is None
                        else bool(np.isclose(actual[k], reference[k], rtol=1e-9, atol=1e-9))
                        for k in reference)
            checks[f"unfiltered_snapshot_match_{iteration}_{species}"] = match
            compared += 1
    checks["independent_snapshot_comparisons_present"] = compared > 0
    report = {"case": metadata["case"], "run": str(run), "region": probe,
              "definition": "T_eV = m/(3*q_e) * weighted mean(|v - weighted_mean(v)|^2)",
              "sample_count_per_species": len(series.iterations), "checks": checks,
              "checks_pass": all(checks.values()), "species": {},
              "analyzer_sha256": digest(Path(__file__)), "inputs_sha256": digest(run / "inputs.txt"),
              "limitations": ["Nonrelativistic second-moment temperature, not a Maxwellian fit",
                  "Finite-region mean flow is removed; unresolved flow gradients can remain",
                  "No Bessel correction: distribution moment, not an unbiased sample estimator",
                  "Particle exchange changes the sampled population; temperature change is not a heat budget",
                  "Beam energy input is tiny and the 10 ns startup does not demonstrate beam heating"]}
    for species in probe["species"]:
        subset = [r for r in rows if r["species"] == species]
        valid = [r["T_eV"] for r in subset if r["T_eV"] is not None]
        report["species"][species] = {"initial_T_eV": subset[0]["T_eV"],
             "final_T_eV": subset[-1]["T_eV"], "min_T_eV": min(valid) if valid else None,
             "max_T_eV": max(valid) if valid else None,
             "min_macro_count": min(r["macro_count"] for r in subset),
             "max_macro_count": max(r["macro_count"] for r in subset)}
    with (run / "temperature_history.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_json(run / "temperature_summary.json", report)
    plot_temperature(rows, metadata, run / "temperature_history.png")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    result = analyze_temperature(args.run_dir.resolve())
    print(json.dumps(result, indent=2))
    if not result["checks_pass"]:
        raise SystemExit(1)
