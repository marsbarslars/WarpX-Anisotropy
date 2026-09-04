#!/usr/bin/env python3
"""Animate the matched 35-degree collisionless/collisional campaign at 24 fps."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import imageio_ffmpeg
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.colors import LogNorm


HERE = Path(__file__).resolve().parent
VIS = HERE / "visualizations"


@dataclass
class CaseData:
    label: str
    color: str
    times_us: np.ndarray
    count: np.ndarray
    survival: np.ndarray
    confined_energy: np.ndarray
    phase: np.ndarray
    xi: np.ndarray
    vparallel_edges: np.ndarray
    vperp_edges: np.ndarray
    xi_edges: np.ndarray
    mirror_ratio: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=VIS / "campaign_velocity_space_035_animation_24fps.mp4",
    )
    parser.add_argument(
        "--preview",
        type=Path,
        default=VIS / "campaign_velocity_space_035_animation_preview.png",
    )
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--duration", type=float, default=12.0)
    parser.add_argument("--dpi", type=int, default=150)
    return parser.parse_args()


def load_case(name: str, label: str, color: str) -> CaseData:
    case = HERE / "runs" / name
    archive = np.load(case / "velocity_space.npz")
    with (case / "timeseries.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    summary = json.loads((case / "campaign_summary.json").read_text(encoding="utf-8"))
    return CaseData(
        label=label,
        color=color,
        times_us=np.array([float(row["time_us"]) for row in rows]),
        count=np.array([float(row["beam_count"]) for row in rows]),
        survival=np.array([float(row["survival_fraction"]) for row in rows]),
        confined_energy=np.array([float(row["confined_energy_fraction"]) for row in rows]),
        phase=archive["vparallel_vperp_counts"].astype(float),
        xi=archive["xi_counts"].astype(float),
        vparallel_edges=archive["vparallel_edges_m_s"].astype(float) / 1000.0,
        vperp_edges=archive["vperp_edges_m_s"].astype(float) / 1000.0,
        xi_edges=archive["xi_edges"].astype(float),
        mirror_ratio=float(summary["field_mirror_ratio"]),
    )


def interpolate(values: np.ndarray, position: float) -> np.ndarray:
    lo = int(np.floor(position))
    hi = min(lo + 1, len(values) - 1)
    fraction = position - lo
    return (1.0 - fraction) * values[lo] + fraction * values[hi]


def probability(values: np.ndarray) -> np.ndarray:
    total = float(np.sum(values))
    return values / total if total > 0.0 else np.zeros_like(values)


def main() -> None:
    args = parse_args()
    if args.fps <= 0 or args.duration <= 0.0:
        raise SystemExit("--fps and --duration must be positive")
    if args.out.suffix.lower() != ".mp4":
        raise SystemExit("--out must end in .mp4")

    control = load_case(
        "campaign_angle_035_collisionless", "collisionless", "#2878b5"
    )
    collision = load_case(
        "campaign_angle_035_collisional", "collisional", "#e26a2c"
    )
    if not np.array_equal(control.times_us, collision.times_us):
        raise SystemExit("The two diagnostics do not share the same time axis")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.preview.parent.mkdir(parents=True, exist_ok=True)
    mpl.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()

    frame_count = int(round(args.fps * args.duration))
    positions = np.linspace(0.0, len(control.times_us) - 1.0, frame_count)
    extent = [
        control.vparallel_edges[0],
        control.vparallel_edges[-1],
        control.vperp_edges[0],
        control.vperp_edges[-1],
    ]
    initial_probability = probability(control.phase[0]).T
    all_phase_probability = np.concatenate(
        [
            control.phase.reshape(len(control.phase), -1)
            / np.maximum(control.phase.sum(axis=(1, 2))[:, None], 1.0),
            collision.phase.reshape(len(collision.phase), -1)
            / np.maximum(collision.phase.sum(axis=(1, 2))[:, None], 1.0),
        ],
        axis=0,
    )
    positive = all_phase_probability[all_phase_probability > 0.0]
    norm = LogNorm(vmin=max(float(positive.min()), 1.0e-5), vmax=float(positive.max()))
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#f7f7f7")

    fig = plt.figure(figsize=(12.8, 7.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(2.15, 1.15))
    phase_axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    pitch_ax = fig.add_subplot(grid[1, 0])
    history_ax = fig.add_subplot(grid[1, 1])

    images = []
    titles = []
    for ax, case in zip(phase_axes, (control, collision)):
        image = ax.imshow(
            np.ma.masked_less_equal(initial_probability, 0.0),
            origin="lower",
            extent=extent,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
        )
        images.append(image)
        titles.append(ax.set_title(""))
        ax.set_xlabel("$v_\\parallel$ [km/s]")
        ax.set_ylabel("$v_\\perp$ [km/s]")
        ax.grid(alpha=0.13)

        x_boundary = np.linspace(0.0, max(abs(extent[0]), abs(extent[1])), 240)
        y_boundary = x_boundary / np.sqrt(case.mirror_ratio - 1.0)
        for sign in (-1.0, 1.0):
            ax.plot(
                sign * x_boundary,
                y_boundary,
                color="#777777",
                linestyle=":",
                linewidth=1.8,
                label="loss-cone boundary" if ax is phase_axes[0] and sign > 0 else None,
            )
        if ax is phase_axes[0]:
            ax.legend(frameon=False, loc="upper left")

    colorbar = fig.colorbar(images[-1], ax=phase_axes, shrink=0.91, pad=0.02)
    colorbar.set_label("Probability per velocity bin")

    xi_center = 0.5 * (control.xi_edges[:-1] + control.xi_edges[1:])
    initial_xi = probability(control.xi[0])
    pitch_ax.step(
        xi_center,
        initial_xi,
        where="mid",
        color="#555555",
        linewidth=2.0,
        label="initial",
    )
    pitch_lines = []
    for case in (control, collision):
        (line,) = pitch_ax.step(
            xi_center,
            initial_xi,
            where="mid",
            color=case.color,
            linewidth=2.1,
            label=case.label,
        )
        pitch_lines.append(line)
    xi_boundary = np.sqrt(1.0 - 1.0 / control.mirror_ratio)
    pitch_ax.axvline(xi_boundary, color="#777777", linestyle=":", linewidth=1.6)
    pitch_ax.axvline(-xi_boundary, color="#777777", linestyle=":", linewidth=1.6)
    pitch_ax.set_xlim(-1.0, 1.0)
    pitch_ax.set_ylim(0.0, max(0.56, 1.08 * initial_xi.max()))
    pitch_ax.set_xlabel("Local pitch cosine $\\xi=v_\\parallel/v$")
    pitch_ax.set_ylabel("Probability")
    pitch_ax.set_title("Pitch redistribution")
    pitch_ax.grid(alpha=0.2)
    pitch_ax.legend(frameon=False, ncol=3, loc="upper left")

    history_lines = []
    for case in (control, collision):
        (survival_line,) = history_ax.plot(
            [], [], color=case.color, linewidth=2.2, label=f"{case.label} survival"
        )
        (energy_line,) = history_ax.plot(
            [],
            [],
            color=case.color,
            linestyle="--",
            linewidth=2.0,
            label=f"{case.label} energy in domain",
        )
        history_lines.append((survival_line, energy_line))
    cursor = history_ax.axvline(0.0, color="#555555", linestyle=":", linewidth=1.5)
    history_ax.set_xlim(control.times_us[0], control.times_us[-1])
    history_ax.set_ylim(0.0, 104.0)
    history_ax.set_xlabel("Time [µs]")
    history_ax.set_ylabel("Fraction of initial beam [%]")
    history_ax.set_title("Survival and beam energy remaining")
    history_ax.grid(alpha=0.2)
    history_ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower left")

    main_title = fig.suptitle("")
    fig.text(
        0.5,
        0.004,
        "24 fps presentation; visual interpolation between 4.4 µs diagnostic snapshots",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#666666",
    )

    def update(frame: int):
        position = positions[frame]
        time_us = float(interpolate(control.times_us, position))
        main_title.set_text(f"35° near-loss-cone evolution · t = {time_us:5.1f} µs")

        artists = [main_title, cursor]
        for image, title, pitch_line, history_pair, case in zip(
            images, titles, pitch_lines, history_lines, (control, collision)
        ):
            phase_probability = probability(interpolate(case.phase, position)).T
            image.set_data(np.ma.masked_less_equal(phase_probability, 0.0))
            current_count = int(round(float(interpolate(case.count, position))))
            current_survival = 100.0 * float(interpolate(case.survival, position))
            current_energy = 100.0 * float(interpolate(case.confined_energy, position))
            mean_survivor_energy = 5.630570733 * current_energy / max(current_survival, 1.0e-12)
            title.set_text(
                f"{case.label.capitalize()}\n"
                f"N = {current_count:,} · survival {current_survival:.1f}% · "
                f"energy in domain {current_energy:.1f}% · mean {mean_survivor_energy:.2f} eV"
            )
            pitch_line.set_ydata(probability(interpolate(case.xi, position)))

            full_index = int(np.floor(position))
            history_time = np.append(case.times_us[: full_index + 1], time_us)
            survival_history = np.append(
                100.0 * case.survival[: full_index + 1], current_survival
            )
            energy_history = np.append(
                100.0 * case.confined_energy[: full_index + 1], current_energy
            )
            history_pair[0].set_data(history_time, survival_history)
            history_pair[1].set_data(history_time, energy_history)
            artists.extend([image, title, pitch_line, *history_pair])
        cursor.set_xdata([time_us, time_us])
        return artists

    preview_frame = int(round(0.68 * (frame_count - 1)))
    update(preview_frame)
    fig.savefig(args.preview, dpi=args.dpi)

    animation = FuncAnimation(
        fig,
        update,
        frames=frame_count,
        interval=1000.0 / args.fps,
        blit=False,
    )
    writer = FFMpegWriter(
        fps=args.fps,
        codec="libx264",
        bitrate=5000,
        extra_args=[
            "-vf",
            "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            "-pix_fmt",
            "yuv420p",
            "-movflags",
            "+faststart",
        ],
    )

    def progress(frame: int, total: int) -> None:
        if frame == 0 or frame + 1 == total or (frame + 1) % args.fps == 0:
            print(f"rendered {frame + 1}/{total} frames")

    animation.save(args.out, writer=writer, dpi=args.dpi, progress_callback=progress)
    plt.close(fig)
    print(f"wrote {args.out}")
    print(f"wrote {args.preview}")


if __name__ == "__main__":
    main()
