#!/usr/bin/env python3
"""Build complete time-series, loss-event, and velocity-space campaign data."""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
TEAM_SCRIPTS = REPO_ROOT / "Plasma-Hackathon" / "scripts"
sys.path.insert(0, str(TEAM_SCRIPTS))

from plasma import load_field_map  # noqa: E402


C = 299792458.0
M_P = 1.67262192595e-27
Q_E = 1.602176634e-19
BOUNDS = np.array([[-1.0, 1.0], [-1.0, 1.0], [0.0, 5.0]])


def snapshot(series, iteration: int, field) -> dict[str, np.ndarray]:
    values = series.get_particle(
        ["x", "y", "z", "id", "ux", "uy", "uz"],
        species="beam_protons",
        iteration=int(iteration),
    )
    x, y, z, ids, ux, uy, uz = (np.asarray(value) for value in values)
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
    velocity = np.column_stack((ux, uy, uz)) * (C / gamma[:, None])
    speed = np.linalg.norm(velocity, axis=1)
    energy_eV = (gamma - 1.0) * M_P * C**2 / Q_E
    position = np.column_stack((x, y, z))
    b = field.interpolate(x, y, z) if len(ids) else np.empty((0, 3))
    bmag = np.linalg.norm(b, axis=1)
    bhat = b / np.where(bmag[:, None] > 0.0, bmag[:, None], np.inf)
    vparallel = np.einsum("ij,ij->i", velocity, bhat)
    vperp = np.sqrt(np.maximum(speed**2 - vparallel**2, 0.0))
    xi = vparallel / np.where(speed > 0.0, speed, np.inf)
    return {
        "id": ids.astype(np.uint64),
        "position": position,
        "velocity": velocity,
        "speed": speed,
        "energy_eV": energy_eV,
        "xi": xi,
        "vparallel": vparallel,
        "vperp": vperp,
    }


def boundary_label(position: np.ndarray, velocity: np.ndarray) -> str:
    names = ("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi")
    distance = np.array(
        [
            position[0] - BOUNDS[0, 0],
            BOUNDS[0, 1] - position[0],
            position[1] - BOUNDS[1, 0],
            BOUNDS[1, 1] - position[1],
            position[2] - BOUNDS[2, 0],
            BOUNDS[2, 1] - position[2],
        ]
    )
    outward_speed = np.array(
        [-velocity[0], velocity[0], -velocity[1], velocity[1], -velocity[2], velocity[2]]
    )
    flight_time = np.where(outward_speed > 0.0, distance / outward_speed, np.inf)
    index = int(np.argmin(flight_time))
    if not np.isfinite(flight_time[index]):
        index = int(np.argmin(distance))
    return names[index]


def first_crossing(rows: list[dict], threshold: float) -> float | None:
    for row in rows:
        if row["survival_fraction"] <= threshold:
            return float(row["time_us"])
    return None


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value):
    """Convert NumPy values and non-finite floats to strict JSON values."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def analyze_case(case: Path, field, loss_cone_deg: float) -> dict:
    from openpmd_viewer import OpenPMDTimeSeries

    manifest = json.loads((case / "manifest.json").read_text(encoding="utf-8"))
    series = OpenPMDTimeSeries(str(case / "diags" / "beam_diag"))
    iterations = np.asarray(series.iterations, dtype=int)
    dt_s = float(manifest["dt_s"])
    xi_loss_boundary = math.cos(math.radians(loss_cone_deg))

    states: list[dict[str, np.ndarray]] = []
    rows: list[dict] = []
    losses: list[dict] = []
    cumulative_lost_energy = 0.0
    previous: dict[str, np.ndarray] | None = None
    initial_count = 0
    initial_energy = 0.0

    for iteration in iterations:
        state = snapshot(series, int(iteration), field)
        states.append(state)
        time_us = iteration * dt_s * 1.0e6
        if previous is None:
            initial_count = len(state["id"])
            initial_energy = float(np.sum(state["energy_eV"]))
        else:
            current_ids = set(state["id"].tolist())
            lost_indices = [i for i, particle_id in enumerate(previous["id"]) if int(particle_id) not in current_ids]
            previous_time_us = rows[-1]["time_us"]
            for index in lost_indices:
                particle_energy = float(previous["energy_eV"][index])
                cumulative_lost_energy += particle_energy
                position = previous["position"][index]
                velocity = previous["velocity"][index]
                losses.append(
                    {
                        "particle_id": int(previous["id"][index]),
                        "loss_time_lower_us": previous_time_us,
                        "loss_time_upper_us": time_us,
                        "loss_time_midpoint_us": 0.5 * (previous_time_us + time_us),
                        "inferred_boundary": boundary_label(position, velocity),
                        "last_x_m": float(position[0]),
                        "last_y_m": float(position[1]),
                        "last_z_m": float(position[2]),
                        "last_energy_eV": particle_energy,
                        "last_xi": float(previous["xi"][index]),
                        "last_vparallel_m_s": float(previous["vparallel"][index]),
                        "last_vperp_m_s": float(previous["vperp"][index]),
                    }
                )

        count = len(state["id"])
        confined_energy = float(np.sum(state["energy_eV"]))
        vparallel2 = float(np.mean(state["vparallel"] ** 2)) if count else math.nan
        vperp2 = float(np.mean(state["vperp"] ** 2)) if count else math.nan
        rows.append(
            {
                "iteration": int(iteration),
                "time_us": time_us,
                "beam_count": count,
                "survival_fraction": count / initial_count if initial_count else math.nan,
                "confined_energy_eV": confined_energy,
                "confined_energy_fraction": confined_energy / initial_energy if initial_energy else math.nan,
                "cumulative_lost_energy_estimate_eV": cumulative_lost_energy,
                "accounted_beam_energy_fraction": (
                    (confined_energy + cumulative_lost_energy) / initial_energy if initial_energy else math.nan
                ),
                "survivor_mean_energy_eV": float(np.mean(state["energy_eV"])) if count else math.nan,
                "survivor_energy_std_eV": float(np.std(state["energy_eV"])) if count else math.nan,
                "mean_xi": float(np.mean(state["xi"])) if count else math.nan,
                "pitch_variance": float(np.var(state["xi"])) if count else math.nan,
                "midplane_loss_cone_fraction": (
                    float(np.mean(np.abs(state["xi"]) >= xi_loss_boundary)) if count else math.nan
                ),
                "velocity_anisotropy_vpar2_over_half_vperp2": (
                    vparallel2 / (0.5 * vperp2) if count and vperp2 > 0.0 else math.nan
                ),
            }
        )
        previous = state

    write_csv(case / "timeseries.csv", rows)
    write_csv(case / "loss_events.csv", losses)

    initial_speed = float(np.median(states[0]["speed"])) if len(states[0]["speed"]) else 1.0
    xi_edges = np.linspace(-1.0, 1.0, 81)
    energy_edges = np.linspace(0.0, max(10.0, 1.5 * float(np.max(states[0]["energy_eV"]))), 81)
    vparallel_edges = np.linspace(-1.5 * initial_speed, 1.5 * initial_speed, 81)
    vperp_edges = np.linspace(0.0, 1.5 * initial_speed, 61)
    xi_counts = np.stack([np.histogram(state["xi"], bins=xi_edges)[0] for state in states])
    energy_counts = np.stack([np.histogram(state["energy_eV"], bins=energy_edges)[0] for state in states])
    phase_counts = np.stack(
        [
            np.histogram2d(
                state["vparallel"], state["vperp"], bins=(vparallel_edges, vperp_edges)
            )[0]
            for state in states
        ]
    )
    np.savez_compressed(
        case / "velocity_space.npz",
        iterations=iterations,
        times_us=np.array([row["time_us"] for row in rows]),
        xi_edges=xi_edges,
        xi_counts=xi_counts,
        energy_edges_eV=energy_edges,
        energy_counts=energy_counts,
        vparallel_edges_m_s=vparallel_edges,
        vperp_edges_m_s=vperp_edges,
        vparallel_vperp_counts=phase_counts,
    )

    reduced_energy = np.loadtxt(
        case / "diags" / "reducedfiles" / "particle_energy.txt",
        comments="#",
        ndmin=2,
    )
    initial_reduced = reduced_energy[0]
    final_reduced = reduced_energy[-1]
    reduced_energy_summary = {
        "initial_beam_total_J": float(initial_reduced[3]),
        "final_beam_total_J": float(final_reduced[3]),
        "initial_background_total_J": float(initial_reduced[4]),
        "final_background_total_J": float(final_reduced[4]),
        "background_energy_change_J": float(final_reduced[4] - initial_reduced[4]),
        "initial_background_mean_eV": float(initial_reduced[7] / Q_E),
        "final_background_mean_eV": float(final_reduced[7] / Q_E),
    }

    final = rows[-1]
    summary = {
        **manifest,
        "field_mirror_ratio": field.mirror_ratio()[2],
        "theoretical_loss_cone_deg": loss_cone_deg,
        "diagnostic_final_time_us": final["time_us"],
        "initial_beam_count": initial_count,
        "final_beam_count": final["beam_count"],
        "final_survival_fraction": final["survival_fraction"],
        "final_confined_energy_fraction": final["confined_energy_fraction"],
        "final_accounted_beam_energy_fraction": final["accounted_beam_energy_fraction"],
        "estimated_collisionally_removed_energy_fraction": 1.0 - final["accounted_beam_energy_fraction"],
        "final_survivor_mean_energy_eV": final["survivor_mean_energy_eV"],
        "final_pitch_variance": final["pitch_variance"],
        "final_midplane_loss_cone_fraction": final["midplane_loss_cone_fraction"],
        "final_velocity_anisotropy": final["velocity_anisotropy_vpar2_over_half_vperp2"],
        "time_to_90pct_survival_us": first_crossing(rows, 0.90),
        "time_to_75pct_survival_us": first_crossing(rows, 0.75),
        "time_to_50pct_survival_us": first_crossing(rows, 0.50),
        "loss_events_recorded": len(losses),
        "losses_by_inferred_boundary": dict(Counter(event["inferred_boundary"] for event in losses)),
        "loss_time_resolution_us": float(np.diff([row["time_us"] for row in rows]).max()) if len(rows) > 1 else None,
        "reduced_energy": reduced_energy_summary,
        "data_products": {
            "time_series": "timeseries.csv",
            "loss_events": "loss_events.csv",
            "velocity_space": "velocity_space.npz",
            "raw_particles": "diags/beam_diag",
        },
    }
    summary = json_safe(summary)
    (case / "campaign_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        f"{case.name}: survival={summary['final_survival_fraction']:.3f}, "
        f"confined E/E0={summary['final_confined_energy_fraction']:.3f}, "
        f"losses={len(losses)}"
    )
    return summary


def main() -> None:
    field = load_field_map(str(HERE / "example-femm-3d.h5"))
    loss_cone_deg = field.loss_cone_deg()
    cases = sorted((HERE / "runs").glob("campaign_angle_*"))
    if not cases:
        raise SystemExit("No campaign cases found; run generate_campaign.py first")
    incomplete = [case for case in cases if not (case / "diags" / "beam_diag").exists()]
    if incomplete:
        raise SystemExit(f"Campaign output missing for: {', '.join(case.name for case in incomplete)}")

    summaries = [analyze_case(case, field, loss_cone_deg) for case in cases]
    compact_keys = [
        "case",
        "angle_deg",
        "collision_model",
        "beam_particles",
        "diagnostic_final_time_us",
        "final_beam_count",
        "final_survival_fraction",
        "final_confined_energy_fraction",
        "final_accounted_beam_energy_fraction",
        "final_survivor_mean_energy_eV",
        "final_pitch_variance",
        "final_midplane_loss_cone_fraction",
        "final_velocity_anisotropy",
        "time_to_90pct_survival_us",
        "time_to_75pct_survival_us",
        "time_to_50pct_survival_us",
        "loss_events_recorded",
    ]
    compact = [{key: summary.get(key) for key in compact_keys} for summary in summaries]
    write_csv(HERE / "campaign_summary.csv", compact)
    (HERE / "campaign_summary.json").write_text(
        json.dumps(json_safe(summaries), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    by_key = {(item["angle_deg"], item["collision_model"]): item for item in summaries}
    comparisons = []
    def difference(left, right):
        return left - right if left is not None and right is not None else None

    for angle in sorted({item["angle_deg"] for item in summaries}):
        control = by_key[(angle, "none")]
        collision = by_key[(angle, "hybrid")]
        comparisons.append(
            {
                "angle_deg": angle,
                "is_above_theoretical_loss_cone": angle > loss_cone_deg,
                "collisionless_survival": control["final_survival_fraction"],
                "collisional_survival": collision["final_survival_fraction"],
                "collision_associated_survival_change": difference(
                    collision["final_survival_fraction"], control["final_survival_fraction"]
                ),
                "collisionless_confined_energy_fraction": control["final_confined_energy_fraction"],
                "collisional_confined_energy_fraction": collision["final_confined_energy_fraction"],
                "collision_associated_confined_energy_change": difference(
                    collision["final_confined_energy_fraction"],
                    control["final_confined_energy_fraction"],
                ),
                "collisionless_pitch_variance": control["final_pitch_variance"],
                "collisional_pitch_variance": collision["final_pitch_variance"],
                "collision_associated_pitch_variance": difference(
                    collision["final_pitch_variance"], control["final_pitch_variance"]
                ),
            }
        )
    write_csv(HERE / "campaign_comparison.csv", comparisons)
    (HERE / "campaign_comparison.json").write_text(
        json.dumps(json_safe(comparisons), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(f"Wrote {HERE / 'campaign_summary.csv'}")
    print(f"Wrote {HERE / 'campaign_comparison.csv'}")


if __name__ == "__main__":
    main()
