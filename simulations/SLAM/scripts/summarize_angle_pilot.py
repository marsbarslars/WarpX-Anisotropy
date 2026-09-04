#!/usr/bin/env python3
"""Combine the moderate-density collisional angle pilot into CSV and JSON."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def main() -> None:
    summary_paths = sorted((HERE / "runs").glob("pilot_angle_*/summary.json"))
    if not summary_paths:
        raise SystemExit("No completed angle-pilot summaries found")

    rows = []
    for path in summary_paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        reservoir = data.get("background_reservoir") or {}
        rows.append(
            {
                "angle_deg": data["angle_deg"],
                "density_m3": data["background_density_m3"],
                "end_time_s": data["end_time_s"],
                "beam_count_initial": data["initial"]["count"],
                "beam_count_final": data["final"]["count"],
                "survival_fraction": data["survival_fraction"],
                "background_retention": reservoir.get("retention_fraction"),
                "energy_initial_eV": data["initial"]["mean_energy_eV"],
                "energy_final_eV": data["final"]["mean_energy_eV"],
                "energy_ratio": data["energy_ratio"],
                "pitch_mean_initial": data["initial"]["mean_xi"],
                "pitch_mean_final": data["final"]["mean_xi"],
                "pitch_std_initial": data["initial"]["std_xi"],
                "pitch_std_final": data["final"]["std_xi"],
            }
        )

    rows.sort(key=lambda row: row["angle_deg"])
    json_path = HERE / "angle_pilot_summary.json"
    csv_path = HERE / "angle_pilot_summary.csv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    for row in rows:
        print(
            f"{row['angle_deg']:>2} deg: survival={row['survival_fraction']:.3f}, "
            f"E/E0={row['energy_ratio']:.6f}, "
            f"<xi>={row['pitch_mean_final']:.6f}, "
            f"sigma_xi={row['pitch_std_final']:.6f}"
        )


if __name__ == "__main__":
    main()
