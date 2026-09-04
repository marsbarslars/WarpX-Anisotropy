#!/usr/bin/env python3
"""Assess whether beam energy has plateaued near the end of a horizon probe."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


Q_E = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--window-points", type=int, default=6)
    parser.add_argument("--plateau-fraction", type=float, default=0.01)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    case = args.case.resolve()
    energy_path = case / "diags" / "reducedfiles" / "particle_energy.txt"
    number_path = case / "diags" / "reducedfiles" / "particle_number.txt"
    energy = np.loadtxt(energy_path, comments="#", ndmin=2)
    number = np.loadtxt(number_path, comments="#", ndmin=2)
    if len(energy) < args.window_points:
        raise SystemExit("Not enough reduced-diagnostic points for the requested window")

    time_s = energy[:, 1]
    mean_energy_eV = energy[:, 6] / Q_E
    window_time = time_s[-args.window_points :]
    window_energy = mean_energy_eV[-args.window_points :]
    slope_eV_s, intercept_eV = np.polyfit(window_time, window_energy, 1)
    window_fractional_change = abs(window_energy[-1] - window_energy[0]) / window_energy[0]
    plateau = bool(window_fractional_change <= args.plateau_fraction)
    survival = number[:, 3] / number[0, 3]

    def first_crossing(threshold: float) -> float | None:
        crossed = np.flatnonzero(survival <= threshold)
        return float(number[crossed[0], 1]) if len(crossed) else None

    report = {
        "case": case.name,
        "time_initial_s": float(time_s[0]),
        "time_final_s": float(time_s[-1]),
        "energy_initial_eV": float(mean_energy_eV[0]),
        "energy_final_eV": float(mean_energy_eV[-1]),
        "energy_ratio": float(mean_energy_eV[-1] / mean_energy_eV[0]),
        "window_points": args.window_points,
        "window_duration_s": float(window_time[-1] - window_time[0]),
        "window_energy_fractional_change": float(window_fractional_change),
        "window_energy_slope_eV_per_s": float(slope_eV_s),
        "linear_intercept_eV": float(intercept_eV),
        "plateau_threshold_fraction": args.plateau_fraction,
        "energy_plateau": plateau,
        "beam_macroparticles_initial": int(number[0, 3]),
        "beam_macroparticles_final": int(number[-1, 3]),
        "background_macroparticles_initial": int(number[0, 4]),
        "background_macroparticles_final": int(number[-1, 4]),
        "survival_fraction_final": float(survival[-1]),
        "time_to_90pct_survival_s": first_crossing(0.90),
        "time_to_75pct_survival_s": first_crossing(0.75),
        "time_to_50pct_survival_s": first_crossing(0.50),
    }
    output = case / "horizon_assessment.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"energy [eV]: {report['energy_initial_eV']:.6f} -> {report['energy_final_eV']:.6f}")
    print(
        f"last {report['window_duration_s'] * 1e6:.1f} us fractional energy change: "
        f"{report['window_energy_fractional_change']:.6f}"
    )
    print(f"energy plateau at {args.plateau_fraction:.1%} threshold: {plateau}")
    print(f"final beam survival: {report['survival_fraction_final']:.3f}")
    if report["time_to_90pct_survival_s"] is not None:
        print(f"90% survival crossed at {report['time_to_90pct_survival_s'] * 1e6:.1f} us")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
