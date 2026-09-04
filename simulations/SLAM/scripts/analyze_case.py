#!/usr/bin/env python3
"""Summarize one hybrid collisional-mirror run from its beam openPMD output."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
TEAM_SCRIPTS = REPO_ROOT / "Plasma-Hackathon" / "scripts"
sys.path.insert(0, str(TEAM_SCRIPTS))

from plasma import load_field_map  # noqa: E402
from plasma.particles import C  # noqa: E402

M_P = 1.67262192595e-27
Q_E = 1.602176634e-19


def load_reduced_bounds(path: Path) -> tuple[np.ndarray, np.ndarray] | None:
    """Return the first and last rows of a WarpX reduced diagnostic."""
    if not path.exists():
        return None
    values = np.loadtxt(path, comments="#", ndmin=2)
    if not len(values):
        return None
    return values[0], values[-1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--field", type=Path, default=HERE / "example-femm-3d.h5")
    return parser.parse_args()


def load_state(series, iteration: int, field) -> dict:
    x, y, z, ids, ux, uy, uz = series.get_particle(
        ["x", "y", "z", "id", "ux", "uy", "uz"],
        species="beam_protons",
        iteration=int(iteration),
    )
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    ux, uy, uz = np.asarray(ux), np.asarray(uy), np.asarray(uz)
    gamma = np.sqrt(1.0 + ux**2 + uy**2 + uz**2)
    velocity = np.column_stack([ux * C / gamma, uy * C / gamma, uz * C / gamma])
    speed = np.linalg.norm(velocity, axis=1)
    energy_eV = 0.5 * M_P * speed**2 / Q_E
    positions = np.column_stack([x, y, z])
    b_local = field.interpolate(x, y, z) if len(x) else np.empty((0, 3))
    bmag = np.linalg.norm(b_local, axis=1)
    xi = (
        np.einsum("ij,ij->i", velocity, b_local)
        / np.where(speed * bmag > 0.0, speed * bmag, np.inf)
        if len(x)
        else np.array([])
    )
    return {
        "count": int(len(x)),
        "ids": np.asarray(ids),
        "mean_energy_eV": float(np.mean(energy_eV)) if len(x) else None,
        "std_energy_eV": float(np.std(energy_eV)) if len(x) else None,
        "mean_xi": float(np.mean(xi)) if len(x) else None,
        "std_xi": float(np.std(xi)) if len(x) else None,
        "mean_position_m": positions.mean(axis=0).tolist() if len(x) else None,
    }


def main() -> None:
    from openpmd_viewer import OpenPMDTimeSeries

    args = parse_args()
    case = args.case.resolve()
    diagnostic = case / "diags" / "beam_diag"
    if not diagnostic.exists():
        raise SystemExit(f"Missing beam diagnostic: {diagnostic}")

    field = load_field_map(str(args.field.resolve()))
    series = OpenPMDTimeSeries(str(diagnostic))
    initial_iteration = int(series.iterations[0])
    final_iteration = int(series.iterations[-1])
    initial = load_state(series, initial_iteration, field)
    final = load_state(series, final_iteration, field)
    initial_ids = initial.pop("ids")
    final_ids = final.pop("ids")
    lost = int(np.sum(~np.isin(initial_ids, final_ids)))

    reduced = case / "diags" / "reducedfiles"
    particle_number = load_reduced_bounds(reduced / "particle_number.txt")
    particle_energy = load_reduced_bounds(reduced / "particle_energy.txt")
    reservoir = None
    reduced_energy = None
    if particle_number is not None:
        number_initial, number_final = particle_number
        reservoir = {
            "initial_macroparticles": int(number_initial[4]),
            "final_macroparticles": int(number_final[4]),
            "retention_fraction": (
                float(number_final[4] / number_initial[4]) if number_initial[4] else None
            ),
            "initial_physical_weight": float(number_initial[7]),
            "final_physical_weight": float(number_final[7]),
        }
    if particle_energy is not None:
        energy_initial, energy_final = particle_energy
        reduced_energy = {
            "initial_beam_mean_eV": float(energy_initial[6] / Q_E),
            "final_beam_mean_eV": float(energy_final[6] / Q_E),
            "initial_background_mean_eV": float(energy_initial[7] / Q_E),
            "final_background_mean_eV": float(energy_final[7] / Q_E),
        }

    manifest = json.loads((case / "manifest.json").read_text(encoding="utf-8"))
    summary = {
        **manifest,
        "initial_iteration": initial_iteration,
        "final_iteration": final_iteration,
        "initial": initial,
        "final": final,
        "lost_particles": lost,
        "survival_fraction": final["count"] / initial["count"] if initial["count"] else None,
        "energy_ratio": (
            final["mean_energy_eV"] / initial["mean_energy_eV"]
            if final["mean_energy_eV"] is not None and initial["mean_energy_eV"]
            else None
        ),
        "pitch_std_change": (
            final["std_xi"] - initial["std_xi"]
            if final["std_xi"] is not None and initial["std_xi"] is not None
            else None
        ),
        "background_reservoir": reservoir,
        "reduced_energy": reduced_energy,
    }
    output = case / "summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"case: {case.name}")
    print(f"beam particles: {initial['count']} -> {final['count']}")
    print(f"lost particles: {lost}")
    print(f"mean beam energy [eV]: {initial['mean_energy_eV']} -> {final['mean_energy_eV']}")
    print(f"pitch std: {initial['std_xi']} -> {final['std_xi']}")
    if reservoir is not None:
        print(
            "background macroparticles: "
            f"{reservoir['initial_macroparticles']} -> {reservoir['final_macroparticles']} "
            f"(retention {reservoir['retention_fraction']:.6f})"
        )
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
