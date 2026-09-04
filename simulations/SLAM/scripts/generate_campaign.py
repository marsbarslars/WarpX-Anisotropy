#!/usr/bin/env python3
"""Generate the matched 300-microsecond mirror angle campaign."""

from __future__ import annotations

from pathlib import Path

from generate_inputs import PARAMETERS, RUNS, write_case


ANGLES_DEG = (15, 25, 30, 35, 40, 50, 60, 75)
DENSITY_M3 = 1.0e16
BEAM_PARTICLES = 5000
RANDOM_SEED = 42
BEAM_DIAG_INTERVAL = 10
BACKGROUND_DIAG_INTERVAL = 50
HORIZON_S = 300.0e-6
# End on a diagnostic step. With dt=0.44 us this is 299.2 us, reported as 300 us.
MAX_STEP = (round(HORIZON_S / PARAMETERS["dt_s"]) // BEAM_DIAG_INTERVAL) * BEAM_DIAG_INTERVAL


def main() -> None:
    cases: list[Path] = []
    for angle in ANGLES_DEG:
        for model in ("none", "hybrid"):
            condition = "collisionless" if model == "none" else "collisional"
            name = f"campaign_angle_{angle:03d}_{condition}"
            cases.append(
                write_case(
                    name,
                    angle_deg=angle,
                    density=DENSITY_M3,
                    collision_model=model,
                    max_step=MAX_STEP,
                    beam_diag_interval=BEAM_DIAG_INTERVAL,
                    background_diag_interval=BACKGROUND_DIAG_INTERVAL,
                    beam_particles=BEAM_PARTICLES,
                    random_seed=RANDOM_SEED,
                )
            )

    launcher = RUNS / "run_campaign.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -e",
        'root="$(cd "$(dirname "$0")" && pwd)"',
    ]
    for case in cases:
        lines.append(f'(cd "$root/{case.name}" && ./run.sh)')
    launcher.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    launcher.chmod(0o755)

    actual_horizon_us = MAX_STEP * PARAMETERS["dt_s"] * 1.0e6
    print(f"Generated {len(cases)} matched campaign cases")
    print(f"Angles: {', '.join(str(a) for a in ANGLES_DEG)} degrees")
    print(f"Horizon: {actual_horizon_us:.2f} microseconds ({MAX_STEP} steps)")
    print(f"Beam macroparticles per case: {BEAM_PARTICLES}")
    print(f"Launcher: {launcher.relative_to(Path(__file__).resolve().parent)}")


if __name__ == "__main__":
    main()
