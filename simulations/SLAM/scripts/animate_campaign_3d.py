#!/usr/bin/env python3
"""Render the matched 35-degree campaign inside the 3D magnetic mirror."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
ANGLE_VIS = REPO_ROOT / "warpx" / "mirror_angle_scan" / "visualization"
sys.path.insert(0, str(ANGLE_VIS))

import animate_angle_scan as base  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=HERE / "visualizations" / "campaign_3d_mirror_035_comparison_24fps.mp4",
    )
    parser.add_argument(
        "--preview",
        type=Path,
        default=HERE / "visualizations" / "campaign_3d_mirror_035_preview.png",
    )
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--duration", type=float, default=12.0)
    parser.add_argument("--n-tracks", type=int, default=120)
    parser.add_argument("--trail", type=int, default=12)
    parser.add_argument("--width", type=int, default=1600)
    parser.add_argument("--height", type=int, default=720)
    return parser.parse_args()


def case_dir(condition: str) -> Path:
    return HERE / "runs" / f"campaign_angle_035_{condition}"


def load_scalar_history(case: Path) -> dict[str, np.ndarray]:
    with (case / "timeseries.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "survival": np.array([float(row["survival_fraction"]) for row in rows]),
        "energy": np.array([float(row["confined_energy_fraction"]) for row in rows]),
    }


def background_points(case: Path) -> np.ndarray:
    series = OpenPMDTimeSeries(str(case / "diags" / "background_sample"))
    iteration = int(series.iterations[0])
    x, y, z = series.get_particle(
        ["x", "y", "z"], species="background_ions", iteration=iteration
    )
    return np.column_stack((x, y, z))


def interpolate_scalar(values: np.ndarray, positions: np.ndarray) -> np.ndarray:
    source = np.arange(len(values), dtype=float)
    return np.interp(positions, source, values)


def expand_tracks(track_set: base.TrackSet, frame_count: int, bounds) -> tuple[base.TrackSet, np.ndarray]:
    """Smooth visible tracks between diagnostics; losses still occur at recorded intervals."""
    source_positions = np.linspace(0.0, len(track_set.iterations) - 1.0, frame_count)
    expanded_position = np.full((frame_count, track_set.positions.shape[1], 3), np.nan)
    expanded_xi = np.full((frame_count, track_set.positions.shape[1]), np.nan)
    xmin, xmax, ymin, ymax, zmin, zmax = bounds

    for output_frame, source_position in enumerate(source_positions):
        lo = int(np.floor(source_position))
        hi = min(lo + 1, len(track_set.iterations) - 1)
        fraction = source_position - lo
        lo_position = track_set.positions[lo]
        hi_position = track_set.positions[hi]
        lo_valid = ~np.isnan(lo_position).any(axis=1)
        hi_valid = ~np.isnan(hi_position).any(axis=1)
        shared = lo_valid & hi_valid
        expanded_position[output_frame, shared] = (
            (1.0 - fraction) * lo_position[shared] + fraction * hi_position[shared]
        )
        expanded_xi[output_frame, shared] = (
            (1.0 - fraction) * track_set.xi[lo, shared]
            + fraction * track_set.xi[hi, shared]
        )

        disappearing = lo_valid & ~hi_valid
        if np.any(disappearing) and fraction < 1.0:
            previous = track_set.positions[max(lo - 1, 0)]
            previous_valid = ~np.isnan(previous).any(axis=1)
            extrapolated = lo_position.copy()
            moving = disappearing & previous_valid
            extrapolated[moving] = lo_position[moving] + fraction * (
                lo_position[moving] - previous[moving]
            )
            inside = (
                (extrapolated[:, 0] >= xmin)
                & (extrapolated[:, 0] <= xmax)
                & (extrapolated[:, 1] >= ymin)
                & (extrapolated[:, 1] <= ymax)
                & (extrapolated[:, 2] >= zmin)
                & (extrapolated[:, 2] <= zmax)
            )
            visible = disappearing & inside
            expanded_position[output_frame, visible] = extrapolated[visible]
            expanded_xi[output_frame, visible] = track_set.xi[lo, visible]

    expanded = base.TrackSet(
        angle_deg=track_set.angle_deg,
        iterations=np.rint(interpolate_scalar(track_set.iterations, source_positions)).astype(int),
        times_s=interpolate_scalar(track_set.times_s, source_positions),
        positions=expanded_position,
        xi=np.clip(expanded_xi, -1.0, 1.0),
        survivor_counts=np.rint(
            interpolate_scalar(track_set.survivor_counts, source_positions)
        ).astype(int),
    )
    return expanded, source_positions


def oblique_camera(plotter, bounds) -> None:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    center = np.array([(xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2])
    span = max(xmax - xmin, ymax - ymin, zmax - zmin)
    camera = center + np.array([0.85 * span, 1.65 * span, 0.42 * span])
    plotter.camera_position = [tuple(camera), tuple(center), (0.0, 0.0, 1.0)]


def render(args: argparse.Namespace) -> None:
    pv = base.import_rendering()
    field = base.load_field_map(str(HERE / "example-femm-3d.h5"))
    grid = base.field_grid(field, pv)
    lines = base.field_lines(grid, pv)
    control_case = case_dir("collisionless")
    collision_case = case_dir("collisional")
    beam_paths = [
        control_case / "diags" / "beam_diag",
        collision_case / "diags" / "beam_diag",
    ]
    keep_ids = base.shared_initial_ids(beam_paths[0], "beam_protons", args.n_tracks)
    raw_tracks = [
        base.load_track_set(path, field, "beam_protons", keep_ids, stride=1)
        for path in beam_paths
    ]
    base.ensure_common_initial_ids(raw_tracks, keep_ids)
    if not np.array_equal(raw_tracks[0].iterations, raw_tracks[1].iterations):
        raise RuntimeError("Campaign diagnostics are not synchronized")

    tracks = raw_tracks
    diagnostic_frames = len(tracks[0].iterations)
    target_video_frames = int(round(args.fps * args.duration))
    if target_video_frames < diagnostic_frames:
        raise ValueError(
            f"Requested duration provides {target_video_frames} video frames for "
            f"{diagnostic_frames} diagnostics; increase --duration"
        )
    repeat_counts = np.full(
        diagnostic_frames, target_video_frames // diagnostic_frames, dtype=int
    )
    extra = target_video_frames - int(repeat_counts.sum())
    if extra:
        repeat_counts[
            np.linspace(0, diagnostic_frames - 1, extra, dtype=int)
        ] += 1
    histories = [load_scalar_history(control_case), load_scalar_history(collision_case)]
    survival = [history["survival"] for history in histories]
    energy = [history["energy"] for history in histories]
    background = background_points(control_case)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.preview.parent.mkdir(parents=True, exist_ok=True)
    plotter = pv.Plotter(
        shape=(1, 2),
        off_screen=True,
        window_size=(args.width, args.height),
        border=True,
    )
    plotter.set_background("#101216")
    conditions = ("Collisionless", "Collisional")
    for column, condition in enumerate(conditions):
        plotter.subplot(0, column)
        prefix = condition.lower()
        base.add_static_scene(
            plotter, grid, lines, pv, prefix, show_bounds=column == 0
        )
        plotter.add_mesh(
            pv.PolyData(background),
            color="#73d2de",
            point_size=3,
            opacity=0.13,
            render_points_as_spheres=False,
            name=f"{prefix}-background",
        )
        oblique_camera(plotter, grid.bounds)
        if column == 0:
            plotter.add_text(
                "cyan: fixed background sample\n"
                "orange: recent beam trajectory\n"
                "sphere color: local pitch ξ (-1 blue, +1 red)",
                position="lower_right",
                font_size=9,
                color="white",
                name="legend",
            )

    plotter.open_movie(str(args.out), framerate=args.fps, quality=9)
    preview_frame = int(round(0.68 * (diagnostic_frames - 1)))
    written_frames = 0
    next_progress = args.fps
    for frame in range(diagnostic_frames):
        for column, (condition, track_set) in enumerate(zip(conditions, tracks)):
            plotter.subplot(0, column)
            prefix = condition.lower()
            trails = base.trail_mesh(track_set, frame, args.trail, pv)
            if trails is None:
                plotter.remove_actor(f"{prefix}-trails", render=False)
            else:
                plotter.add_mesh(
                    trails,
                    color="#ffb000",
                    line_width=1.6,
                    opacity=0.78,
                    name=f"{prefix}-trails",
                )
            heads = base.head_mesh(track_set, frame, pv)
            if heads is None:
                plotter.remove_actor(f"{prefix}-heads", render=False)
            else:
                plotter.add_mesh(
                    heads,
                    scalars="xi",
                    cmap="coolwarm",
                    clim=(-1.0, 1.0),
                    point_size=10,
                    render_points_as_spheres=True,
                    show_scalar_bar=False,
                    name=f"{prefix}-heads",
                )
            plotter.add_text(
                f"{condition} · 35° injection\n"
                f"time: {track_set.times_s[frame] * 1.0e6:5.1f} µs\n"
                f"beam survival: {100.0 * survival[column][frame]:5.1f}% "
                f"({track_set.survivor_counts[frame]:,}/5,000)\n"
                f"confined beam energy: {100.0 * energy[column][frame]:5.1f}%",
                position="upper_left",
                font_size=11,
                color="white",
                name=f"{prefix}-status",
            )
        if frame == preview_frame:
            plotter.screenshot(str(args.preview))
        for _ in range(int(repeat_counts[frame])):
            plotter.write_frame()
            written_frames += 1
            if written_frames >= next_progress:
                print(f"rendered {written_frames}/{target_video_frames} frames")
                next_progress += args.fps
    plotter.close()
    print(f"wrote {args.out}")
    print(f"wrote {args.preview}")


def main() -> None:
    args = parse_args()
    if args.fps <= 0 or args.duration <= 0.0 or args.n_tracks <= 0:
        raise SystemExit("fps, duration, and n-tracks must be positive")
    if args.out.suffix.lower() != ".mp4" or args.preview.suffix.lower() != ".png":
        raise SystemExit("--out must be .mp4 and --preview must be .png")
    render(args)


if __name__ == "__main__":
    main()
