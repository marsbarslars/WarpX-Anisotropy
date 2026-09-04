"""Loading WarpX/openPMD field maps.

The one thing worth knowing: WarpX writes mesh arrays in *file* axis order,
which for these runs is (z, y, x) — not (x, y, z). Reading the array without
consulting ``axisLabels`` silently profiles the wrong axis, which looks like a
physics error rather than an indexing one. Everything here normalises to
(x, y, z) on load.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class FieldMap:
    """A 3D vector field on a uniform grid, axes ordered (x, y, z)."""

    bx: np.ndarray
    by: np.ndarray
    bz: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    @property
    def magnitude(self) -> np.ndarray:
        return np.sqrt(self.bx**2 + self.by**2 + self.bz**2)

    def on_axis_bz(self) -> tuple[np.ndarray, np.ndarray]:
        """(z, Bz) along the machine axis, nearest grid line to x=y=0."""
        i = int(np.argmin(np.abs(self.x)))
        j = int(np.argmin(np.abs(self.y)))
        return self.z, self.bz[i, j, :]

    def mirror_ratio(self) -> tuple[float, float, float]:
        """(B_min, B_max, ratio) on axis — what sets the loss cone."""
        _, bz = self.on_axis_bz()
        lo, hi = float(np.abs(bz).min()), float(np.abs(bz).max())
        return lo, hi, hi / lo

    def loss_cone_deg(self) -> float:
        """Half-angle of the loss cone for particles born at the midplane."""
        lo, hi, _ = self.mirror_ratio()
        return float(np.degrees(np.arcsin(np.sqrt(lo / hi))))

    def interpolate(self, px, py, pz) -> np.ndarray:
        """Trilinear B at particle positions. Returns (N, 3)."""
        out = np.empty((len(px), 3))
        for k, comp in enumerate((self.bx, self.by, self.bz)):
            out[:, k] = _trilinear(comp, self.x, self.y, self.z, px, py, pz)
        return out


def _trilinear(vals, x, y, z, px, py, pz):
    def axis(a, p):
        dx = a[1] - a[0]
        f = np.clip((p - a[0]) / dx, 0, len(a) - 1.000001)
        i = f.astype(int)
        return i, f - i

    ix, fx = axis(x, px)
    iy, fy = axis(y, py)
    iz, fz = axis(z, pz)

    out = np.zeros(len(px))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((fx if dx else 1 - fx)
                     * (fy if dy else 1 - fy)
                     * (fz if dz else 1 - fz))
                out += w * vals[ix + dx, iy + dy, iz + dz]
    return out


def load_field_map(series_path: str, iteration: int | None = None) -> FieldMap:
    """Read the B field of one iteration from a WarpX openPMD diagnostic."""
    from openpmd_viewer import OpenPMDTimeSeries

    ts = OpenPMDTimeSeries(series_path)
    if iteration is None:
        iteration = ts.iterations[0]

    comps, info = {}, None
    for c in ("x", "y", "z"):
        comps[c], info = ts.get_field(field="B", coord=c, iteration=iteration)

    # openpmd-viewer reports the axis order it found; normalise to (x, y, z).
    axes = [info.axes[i] for i in range(3)]
    order = [axes.index(a) for a in ("x", "y", "z")]
    return FieldMap(
        bx=np.transpose(comps["x"], order),
        by=np.transpose(comps["y"], order),
        bz=np.transpose(comps["z"], order),
        x=info.x, y=info.y, z=info.z,
    )
