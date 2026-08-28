#!/usr/bin/env python3
"""
Animate proton trajectories in the WarpX magnetic mirror example.

Run from the simulation directory (the one containing diags/).

    python animate_mirror.py                 # interactive window
    python animate_mirror.py --gif out.gif   # write an animated GIF
    python animate_mirror.py --mp4 out.mp4   # needs: uv pip install imageio-ffmpeg

Magnetic field lines are drawn by default, seeded on a mid-plane ring and
traced both ways; pass --isosurfaces for |B| contours instead.

The mirror axis (z) is rendered horizontally.
"""

import argparse
import os
import sys

import numpy as np
import pyvista as pv
from openpmd_viewer import OpenPMDTimeSeries

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plasma import load_field_map, load_particles          # noqa: E402
from plasma.charts import ChartPanel, stack_side_by_side   # noqa: E402

# ---------------------------------------------------------------- parameters

DIAG_PATH = "diags/diag1"
N_TRACKS = 60      # particles to draw; all 1000 is unreadable
STRIDE = 2         # use every Nth iteration
TRAIL = 80         # trail length in frames; None = full history
ISOSURFACES = 4    # |B| contour levels


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gif", metavar="PATH")
    p.add_argument("--mp4", metavar="PATH")
    p.add_argument("--path", default=DIAG_PATH)
    p.add_argument("--n-tracks", type=int, default=N_TRACKS)
    p.add_argument("--stride", type=int, default=STRIDE)
    p.add_argument("--species", default=None,
                   help="species name (default: the only one in the series)")
    p.add_argument("--png", metavar="PATH",
                   help="render a single still of the final frame")
    p.add_argument("--trail", type=int, default=TRAIL,
                   help="trail length in frames; 0 for full history")
    p.add_argument("--zoom", type=float, default=1.0,
                   help=">1 tightens the framing around the domain")
    p.add_argument("--opacity", type=float, default=0.55,
                   help="opacity of the field rendering")
    p.add_argument("--isosurfaces", action="store_true",
                   help="draw |B| isosurfaces instead of field lines")
    p.add_argument("--n-lines", type=int, default=48,
                   help="approximate number of field lines")
    p.add_argument("--line-width", type=float, default=1.6)
    p.add_argument("--no-charts", action="store_true",
                   help="render only the 3D view, without the live panel")
    p.add_argument("--chart-width", type=int, default=560)
    return p.parse_args()


def resolve_species(ts, requested):
    """Species names differ between runs (protons, ions, ...); pick sensibly."""
    available = list(ts.avail_species or [])
    if requested:
        if requested not in available:
            raise SystemExit(
                f"species {requested!r} not in {available}")
        return requested
    if len(available) == 1:
        return available[0]
    raise SystemExit(f"--species required; found {available}")


# ------------------------------------------------------------------ B field

def load_field(ts, iteration):
    """PyVista ImageData carrying both |B| and the B vector, axes (x, y, z).

    The vector is what field-line tracing needs; |B| colours it.
    """
    comps = {}
    for c in ("x", "y", "z"):
        comps[c], info = ts.get_field(field="B", coord=c, iteration=iteration)

    # openPMD stores array axes in file order, not necessarily (x, y, z).
    axes = [info.axes[i] for i in range(3)]
    order = [axes.index(a) for a in ("x", "y", "z")]
    bx, by, bz = (np.transpose(comps[c], order) for c in ("x", "y", "z"))
    bmag = np.sqrt(bx**2 + by**2 + bz**2)

    coords = {a: getattr(info, a) for a in ("x", "y", "z")}
    origin = tuple(coords[a][0] for a in ("x", "y", "z"))
    spacing = tuple(
        (coords[a][1] - coords[a][0]) if len(coords[a]) > 1 else 1.0
        for a in ("x", "y", "z")
    )

    grid = pv.ImageData(dimensions=bmag.shape, origin=origin, spacing=spacing)
    # ImageData point order is x-fastest, so every array is ravelled Fortran-wise.
    grid["Bmag"] = bmag.ravel(order="F")
    grid["B"] = np.column_stack([c.ravel(order="F") for c in (bx, by, bz)])
    return grid


def field_lines(grid, n_lines, seed_radius_frac=0.75):
    """Trace B from a ring of seeds on the mid-plane, in both directions.

    Seeding on a mid-plane ring rather than a volume is what makes a mirror
    legible: every line runs the length of the machine and converges at the
    throats, instead of a thicket of short segments.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = grid.bounds
    cx, cy, cz = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2
    r = seed_radius_frac * min(xmax - xmin, ymax - ymin) / 2

    # A few concentric rings, so the flux surfaces nest visibly.
    pts = []
    for frac in (0.25, 0.55, 1.0):
        n = max(3, int(n_lines * frac / 1.8))
        theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
        pts.append(np.column_stack([
            cx + frac * r * np.cos(theta),
            cy + frac * r * np.sin(theta),
            np.full(n, cz),
        ]))
    seeds = pv.PolyData(np.vstack(pts))

    span = max(xmax - xmin, ymax - ymin, zmax - zmin)
    return grid.streamlines_from_source(
        seeds, vectors="B",
        integration_direction="both",
        max_length=6 * span,
        initial_step_length=0.2,
        terminal_speed=1e-12,
    )


# --------------------------------------------------------------- trajectories

def load_tracks(ts, iterations, n_tracks, species):
    """(n_frames, n_tracks, 3) array of positions; NaN where absorbed."""
    ids0 = ts.get_particle(["id"], species=species, iteration=iterations[0])[0]
    keep = np.sort(ids0)[:n_tracks]          # stable subset across frames
    order = np.argsort(keep)
    sorted_keep = keep[order]

    tracks = np.full((len(iterations), len(keep), 3), np.nan)

    for k, it in enumerate(iterations):
        x, y, z, ids = ts.get_particle(
            ["x", "y", "z", "id"], species=species, iteration=it
        )
        pos = np.searchsorted(sorted_keep, ids)
        pos_c = np.clip(pos, 0, len(sorted_keep) - 1)
        hit = sorted_keep[pos_c] == ids
        tracks[k, order[pos_c[hit]]] = np.column_stack([x[hit], y[hit], z[hit]])

    return tracks


def trail_mesh(tracks, frame, trail):
    """Polyline mesh of each particle's recent history up to `frame`."""
    start = 0 if trail is None else max(0, frame - trail)
    blocks = []
    for t in tracks[start : frame + 1].transpose(1, 0, 2):
        good = t[~np.isnan(t).any(axis=1)]
        if len(good) > 1:
            blocks.append(pv.lines_from_points(good))
    if not blocks:
        return None
    return pv.MultiBlock(blocks).combine()


def head_mesh(tracks, frame):
    pts = tracks[frame]
    pts = pts[~np.isnan(pts).any(axis=1)]
    return pv.PolyData(pts) if len(pts) else None


# ------------------------------------------------------------------- camera

def side_on_camera(plotter, bounds, zoom=1.0, window=(1400, 700)):
    """Look down -y so the mirror axis (z) runs horizontally across the view.

    Parallel projection, sized to the domain: a perspective view of a long thin
    machine draws the near and far faces of the bounding box at visibly
    different sizes, which reads as field lines escaping the domain when they
    are simply closer to the camera.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    cx, cy, cz = (xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2
    span = max(xmax - xmin, zmax - zmin)
    plotter.camera_position = [
        (cx, cy - 2.5 * span, cz),   # eye
        (cx, cy, cz),                # focal point
        (1.0, 0.0, 0.0),             # view up = +x  ->  z is horizontal
    ]
    plotter.enable_parallel_projection()

    # Vertical on screen is x, horizontal is z. Fit whichever needs more room.
    aspect = window[0] / window[1]
    half = max((xmax - xmin) / 2, (zmax - zmin) / 2 / aspect)
    plotter.camera.parallel_scale = half * 1.06 / zoom


# ---------------------------------------------------------------- compositing

def composite(plotter, panel, fmap, path, species, iteration, ts, n_total):
    """One finished frame: 3D render on the left, live charts on the right."""
    scene = plotter.screenshot(return_img=True)

    part = load_particles(path, iteration, species, ts=ts)
    if len(part):
        b = fmap.interpolate(part.x, part.y, part.z)
        vpar, vperp = part.v_par_perp(b)
        bmag = np.linalg.norm(b, axis=1)
    else:
        vpar = vperp = bmag = np.array([])

    return stack_side_by_side(scene, panel.render(bmag, vpar, vperp, n_total))


# --------------------------------------------------------------------- main

def main():
    args = parse_args()

    ts = OpenPMDTimeSeries(args.path)
    iterations = list(ts.iterations)[:: args.stride]
    if len(iterations) < 2:
        raise SystemExit(
            f"Only {len(iterations)} iteration(s) in {args.path}. "
            "Re-run WarpX with diag1.intervals = 1."
        )

    species = resolve_species(ts, args.species)
    print(f"{len(iterations)} frames, {args.n_tracks} tracks, species {species!r}")
    grid = load_field(ts, iterations[0])
    tracks = load_tracks(ts, iterations, args.n_tracks, species)

    off_screen = bool(args.gif or args.mp4 or args.png)
    charts = not args.no_charts and off_screen
    # 768 is divisible by 16, which keeps ffmpeg from silently resizing.
    win = (1200, 768) if charts else (1400, 700)
    p = pv.Plotter(off_screen=off_screen, window_size=win)
    p.set_background("#07090F")   # matches plasma.charts.GROUND

    if args.isosurfaces:
        p.add_mesh(
            grid.contour(isosurfaces=ISOSURFACES, scalars="Bmag"),
            scalars="Bmag", cmap="cividis", opacity=args.opacity,
            show_scalar_bar=False, name="field",
        )
    else:
        lines = field_lines(grid, args.n_lines)
        if lines.n_points:
            p.add_mesh(
                lines, scalars="Bmag", cmap="cividis",
                opacity=args.opacity, line_width=args.line_width,
                show_scalar_bar=False, name="field",
            )
        else:
            print("warning: no field lines traced; falling back to isosurfaces",
                  file=__import__("sys").stderr)
            p.add_mesh(grid.contour(isosurfaces=ISOSURFACES, scalars="Bmag"),
                       scalars="Bmag", cmap="cividis", opacity=0.3,
                       show_scalar_bar=False, name="field")
    p.add_mesh(pv.Box(grid.bounds), style="wireframe",
               color="gray", opacity=0.3, name="box")

    side_on_camera(p, grid.bounds, args.zoom, win)

    panel = None
    if charts:
        fmap = load_field_map(args.path, iterations[0])
        b_lo, b_hi, _ = fmap.mirror_ratio()
        p0 = load_particles(args.path, iterations[0], species, ts=ts)
        v_scale = float(np.percentile(p0.speed, 99.5) / 1e6)
        panel = ChartPanel(args.chart_width, win[1], (b_lo, b_hi), v_scale,
                           fmap.loss_cone_deg())
        n_total = len(p0)

        # Calibrate the density scale on frames spread across the run, not on
        # frame 0 alone: the distribution concentrates as the loss cone empties.
        probe = iterations[:: max(1, len(iterations) // 6)][:6]
        samples = []
        for it in probe:
            q = load_particles(args.path, it, species, ts=ts)
            if len(q):
                bq = fmap.interpolate(q.x, q.y, q.z)
                samples.append(q.v_par_perp(bq))
        panel.calibrate(samples)

        print(f"charts: |B| {b_lo:.3f}-{b_hi:.3f} T, loss cone "
              f"{fmap.loss_cone_deg():.1f} deg, {n_total} particles, "
              f"density ceiling calibrated on {len(samples)} frames")

    writer = None
    if charts and (args.gif or args.mp4):
        import imageio.v2 as imageio
        if args.gif:
            writer = imageio.get_writer(args.gif, fps=25)
        else:
            writer = imageio.get_writer(args.mp4, fps=30,
                                        macro_block_size=None)

    if charts:
        pass                          # frames are composited and written below
    elif args.gif:
        p.open_gif(args.gif, fps=25)
    elif args.mp4:
        p.open_movie(args.mp4, framerate=30)
    elif args.png:
        pass
    else:
        p.show(interactive_update=True, auto_close=False)

    for frame in range(len(tracks)):
        trails = trail_mesh(tracks, frame, args.trail or None)
        if trails is not None:
            p.add_mesh(trails, color="orange", line_width=1.5,
                       opacity=0.8, name="trails")

        heads = head_mesh(tracks, frame)
        if heads is not None:
            p.add_mesh(heads, color="red", point_size=8,
                       render_points_as_spheres=True, name="heads")

        p.add_text(f"step {iterations[frame]}", position="upper_left",
                   font_size=10, color="white", name="label")

        if charts:
            frame_img = composite(p, panel, fmap, args.path, species,
                                  iterations[frame], ts, n_total)
            if writer is not None:
                writer.append_data(frame_img)
            last_frame = frame_img
        elif args.png:
            pass                      # only the final frame is kept
        elif off_screen:
            p.write_frame()
        else:
            p.update()

    if charts:
        if writer is not None:
            writer.close()
            print(f"wrote {args.gif or args.mp4}")
        if args.png:
            import imageio.v2 as imageio
            imageio.imwrite(args.png, last_frame)
            print(f"wrote {args.png}")
        panel.close()
        p.close()
    elif args.png:
        p.screenshot(args.png)
        p.close()
        print(f"wrote {args.png}")
    elif off_screen:
        p.close()
        print(f"wrote {args.gif or args.mp4}")
    else:
        p.show()


if __name__ == "__main__":
    main()
