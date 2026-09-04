#!/usr/bin/env python3
"""Compare each collisional angle pilot with its matched mirror-only control."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"


def read_summary(pattern: str) -> dict[int, dict]:
    found = {}
    for path in RUNS.glob(pattern):
        summary_path = path / "summary.json"
        if summary_path.exists():
            data = json.loads(summary_path.read_text(encoding="utf-8"))
            found[int(data["angle_deg"])] = data
    return found


def main() -> None:
    hybrid = read_summary("pilot_angle_*_n_1e16_hybrid")
    control = read_summary("control_angle_*_n_1e16_none")
    common = sorted(set(hybrid) & set(control))
    if not common:
        raise SystemExit("No matched hybrid/control summaries found")

    rows = []
    for angle in common:
        coll = hybrid[angle]
        base = control[angle]
        coll_std = coll["final"]["std_xi"]
        base_std = base["final"]["std_xi"]
        rows.append(
            {
                "angle_deg": angle,
                "control_survival": base["survival_fraction"],
                "hybrid_survival": coll["survival_fraction"],
                "collision_associated_loss_fraction": (
                    base["survival_fraction"] - coll["survival_fraction"]
                ),
                "control_energy_ratio": base["energy_ratio"],
                "hybrid_energy_ratio": coll["energy_ratio"],
                "collision_associated_energy_loss_fraction": (
                    base["energy_ratio"] - coll["energy_ratio"]
                ),
                "control_pitch_mean_final": base["final"]["mean_xi"],
                "hybrid_pitch_mean_final": coll["final"]["mean_xi"],
                "collision_associated_pitch_mean_shift": (
                    coll["final"]["mean_xi"] - base["final"]["mean_xi"]
                ),
                "control_pitch_std_final": base_std,
                "hybrid_pitch_std_final": coll_std,
                "collision_associated_pitch_variance": coll_std**2 - base_std**2,
            }
        )

    json_path = HERE / "angle_pilot_comparison.json"
    csv_path = HERE / "angle_pilot_comparison.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    for row in rows:
        print(
            f"{row['angle_deg']:>2} deg: collision-associated d(E/E0)="
            f"{-row['collision_associated_energy_loss_fraction']:.6f}, "
            f"d<xi>={row['collision_associated_pitch_mean_shift']:.6f}, "
            f"excess pitch variance={row['collision_associated_pitch_variance']:.6f}"
        )


if __name__ == "__main__":
    main()
