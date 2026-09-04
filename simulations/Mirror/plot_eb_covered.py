#!/usr/bin/env python3
"""Plot the cell-centered eb_covered field on an XY plane."""

import argparse
import os
import sys

import numpy as np
import pyvista as pv
from openpmd_viewer import OpenPMDTimeSeries


DEFAULT_PATH = "diags/diag1"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", default=DEFAULT_PATH,
                        help="OpenPMD diagnostic path")
    parser.add_argument("--iteration", type=int, default=None,
                        help="iteration to plot (default: first available)")
    parser.add_argument("--z-index", type=int, default=None,
                        help="Z cell index (default: middle cell layer)")
    parser.add_argument("--output", metavar="PATH",
                        help="save a screenshot instead of opening a window")
    return parser.parse_args()


def load_xy_slice(path, iteration, z_index):
    ts = OpenPMDTimeSeries(path)
    available = list(ts.iterations)
    if not available:
        raise SystemExit(f"No iterations found in {path}")
    if iteration is None:
        iteration = int(available[0])
    elif iteration not in available:
        raise SystemExit(f"iteration {iteration} not in {available}")

    values, info = ts.get_field(field="eb_covered", iteration=iteration)

    print("eb_covered min:", np.min(values))
    print("eb_covered max:", np.max(values))

    unique = np.unique(values)
    print("Number of unique values:", len(unique))
    print("First 20 unique values:", unique[:20])

    axes = [info.axes[i] for i in range(values.ndim)]
    order = [axes.index(axis) for axis in ("x", "y", "z")]
    values = np.transpose(values, order)

    nx, ny, nz = values.shape
    if z_index is None:
        z_index = nz // 2
    if not 0 <= z_index < nz:
        raise SystemExit(f"--z-index must be between 0 and {nz - 1}")

    coords = {axis: np.asarray(getattr(info, axis))
              for axis in ("x", "y", "z")}
    spacing = tuple(
        float(coords[axis][1] - coords[axis][0])
        if len(coords[axis]) > 1 else 1.0
        for axis in ("x", "y", "z")
    )
    origin = tuple(float(coords[axis][0] - spacing[i] / 2)
                   for i, axis in enumerate(("x", "y", "z")))

    grid = pv.ImageData(
        dimensions=(nx + 1, ny + 1, 2),
        origin=(origin[0], origin[1], origin[2] + z_index * spacing[2]),
        spacing=(spacing[0], spacing[1], spacing[2]),
    )
    grid.cell_data["eb_covered"] = values[:, :, z_index].ravel(order="F")
    return grid, iteration, z_index


def main():
    args = parse_args()
    path = args.path
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)

    grid, iteration, z_index = load_xy_slice(path, args.iteration, args.z_index)
    plotter = pv.Plotter(off_screen=args.output is not None, window_size=(1000, 800))
    plotter.set_background("#07090F")
    plotter.add_mesh(
        grid,
        scalars="eb_covered",
        cmap="viridis",
        clim=(0, 1),
        show_edges=True,
        edge_color="#526070",
        line_width=0.5,
        show_scalar_bar=True,
        scalar_bar_args={
            "title": "eb_covered",
            "n_labels": 6,
        },
    )
    plotter.view_xy()
    plotter.add_text(
        f"eb_covered, iteration {iteration}, z cell {z_index}",
        position="upper_left", color="white", font_size=12,
    )

    if args.output:
        plotter.screenshot(args.output)
        print(f"wrote {args.output}")
    else:
        plotter.show()


if __name__ == "__main__":
    main()
