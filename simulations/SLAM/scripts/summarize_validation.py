#!/usr/bin/env python3
"""Combine completed collisional-mirror validation summaries into CSV and JSON."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"


def main() -> None:
    summary_paths = sorted(RUNS.glob("validation_*/summary.json"))
    smoke = RUNS / "smoke_fixed_angle_035_n_1e17" / "summary.json"
    if smoke.exists():
        summary_paths.append(smoke)
    if not summary_paths:
        raise SystemExit("No completed validation summaries found")

    rows = []
    for path in summary_paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        reservoir = data.get("background_reservoir") or {}
        rows.append(
            {
                "case": data["case"],
                "density_m3": data["background_density_m3"],
                "collision_model": data.get("collision_model", "hybrid"),
                "dt_s": data["dt_s"],
                "end_time_s": data["end_time_s"],
                "beam_count_initial": data["initial"]["count"],
                "beam_count_final": data["final"]["count"],
                "background_retention": reservoir.get("retention_fraction"),
                "energy_initial_eV": data["initial"]["mean_energy_eV"],
                "energy_final_eV": data["final"]["mean_energy_eV"],
                "energy_ratio": data["energy_ratio"],
                "pitch_std_initial": data["initial"]["std_xi"],
                "pitch_std_final": data["final"]["std_xi"],
                "pitch_std_change": data.get("pitch_std_change"),
            }
        )

    json_path = HERE / "validation_summary.json"
    csv_path = HERE / "validation_summary.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    for row in rows:
        print(
            f"{row['case']}: E/E0={row['energy_ratio']:.6f}, "
            f"sigma_xi={row['pitch_std_final']:.6f}, "
            f"background_retention={row['background_retention']:.6f}"
        )


if __name__ == "__main__":
    main()
