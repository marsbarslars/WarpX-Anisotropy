#!/usr/bin/env python3
"""Beam-energy scan for a SLAM racetrack straight leg.

From inside runs/racetrack-leg:

    ../../.venv/bin/python ../../scripts/racetrack_scan.py --plot racetrack.png

The device sets a hard geometric limit that has nothing to do with the loss
cone: a fast ion's gyroradius has to fit the vessel bore. At the leg midplane
B = 0.103 T, so a 30 keV deuteron gyrates with r = 0.34 m in a 0.23 m bore and
walks into the wall regardless of its pitch.
"""

from __future__ import annotations

import argparse
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plasma import confined_fraction, load_scraped  # noqa: E402

E_KEV = [3, 8, 15, 30]
Q = 1.602176634e-19
M_D = 2.014 * 1.66053907e-27
BORE = 0.230          # m, measured from SLAM_VV.stl
FIELD_FILE = "SLAM_vC5_warpX.h5"
STL_FILE = "SLAM_VV.stl"


def closed_axis_field(path=FIELD_FILE, seed=(0.0, 0.530, 0.0), ds=0.004):
    """Trace the magnetic axis around the loop; return (arclength, |B|).

    A straight line through the machine is not the axis: outside the straight
    sections it leaves the plasma and passes near coil filaments, which turns
    the profile into spikes and makes the mirror ratio meaningless. Following
    the field line is the only way to get |B| along the orbit particles see.
    """
    import h5py
    f = h5py.File(path, "r")
    B = f["data/0/meshes/B"]
    off, sp = B.attrs["gridGlobalOffset"], B.attrs["gridSpacing"]
    comp = [B[k][:] for k in "xyz"]
    n = comp[0].shape
    ax = [off[i] + sp[i] * np.arange(n[i]) for i in range(3)]

    def at(p):
        idx, fr = [], []
        for a, g in zip(p, ax):
            t = (a - g[0]) / (g[1] - g[0])
            if t < 0 or t > len(g) - 1:
                return None
            t = min(t, len(g) - 1.001)
            i = int(t); idx.append(i); fr.append(t - i)
        out = []
        for arr in comp:
            v = 0.0
            for dx in (0, 1):
                for dy in (0, 1):
                    for dz in (0, 1):
                        w = ((fr[0] if dx else 1 - fr[0])
                             * (fr[1] if dy else 1 - fr[1])
                             * (fr[2] if dz else 1 - fr[2]))
                        v += w * arr[idx[0] + dx, idx[1] + dy, idx[2] + dz]
            out.append(v)
        return np.array(out)

    p0 = np.array(seed, float); p = p0.copy()
    S, Bs, D, L = [0.0], [np.linalg.norm(at(p))], [0.0], 0.0
    for _ in range(4000):
        b = at(p)
        if b is None:
            break
        m = np.linalg.norm(b)
        if m < 1e-9:
            break
        p = p + ds * b / m; L += ds
        S.append(L); Bs.append(m); D.append(float(np.linalg.norm(p - p0)))
    S, Bs, D = np.array(S), np.array(Bs), np.array(D)

    # Closure: the first return to the seed after leaving its neighbourhood.
    # A fixed tolerance misses it, because the trace steps past the seed by up
    # to `ds` and need not land inside any given radius.
    away = S > 1.0
    if away.any():
        i = int(np.arange(len(S))[away][np.argmin(D[away])])
        return S[:i + 1], Bs[:i + 1]
    return S, Bs


def vessel_outline(path=STL_FILE):
    with open(path, "rb") as fh:
        fh.read(80)
        nb = struct.unpack("<I", fh.read(4))[0]
        d = np.frombuffer(fh.read(),
                          dtype=np.dtype([("n", "<3f4"), ("v", "<3,3f4"),
                                          ("a", "<u2")]), count=nb)
    return d["v"].astype(np.float64) / 1000.0     # mm -> m


def analyse(energies):
    rows = []
    for E in energies:
        t = f"E{E:g}"
        if not os.path.isdir(f"diags/{t}_reduced"):
            continue
        _, frac = confined_fraction(f"diags/{t}_reduced")
        s, w = load_scraped(f"diags/{t}_scraped", "ions", ())
        cnt = (lambda keys: int(np.isin(w, keys).sum())) if s is not None else (lambda k: 0)
        v = np.sqrt(2 * E * 1e3 * Q / M_D)
        rows.append(dict(E=E, final=float(frac[-1]) if frac is not None else np.nan,
                         wall=cnt(["eb"]),
                         bend=cnt(["xlo", "xhi"]),
                         zedge=cnt(["zlo", "zhi", "ylo", "yhi"]),
                         rg=M_D * v / (Q * 0.1033)))
    return rows


def make_figure(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    E = np.array([r["E"] for r in rows], float)
    final = np.array([r["final"] for r in rows])
    wall = np.array([r["wall"] for r in rows], float)
    bend = np.array([r["bend"] for r in rows], float)
    zedge = np.array([r["zedge"] for r in rows], float)
    n0 = wall + bend + zedge + final * 20000

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9))
    fig.suptitle("SLAM racetrack · full closed loop, real field and vessel",
                 fontsize=13, y=0.98)

    # -- A: the vessel we actually loaded --------------------------------
    ax = axes[0, 0]
    tris = vessel_outline()
    segs = []
    for t in tris[::3]:
        p = t[:, :2]
        segs += [[p[0], p[1]], [p[1], p[2]], [p[2], p[0]]]
    ax.add_collection(LineCollection(np.array(segs), linewidths=0.12,
                                     colors="#8894A8", alpha=0.6))
    ax.add_patch(plt.Rectangle((-1.40, -0.80), 2.80, 1.60, fill=False,
                               edgecolor="#1F6F4A", lw=1.8, ls="--"))
    ax.text(0.0, 0.0, "modelled region\n(whole loop)", ha="center", va="center",
            color="#1F6F4A", fontsize=9)
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-0.95, 0.95)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title("A · The whole racetrack is now modelled", loc="left", fontsize=11)

    # -- B: field along the leg ------------------------------------------
    ax = axes[0, 1]
    S, bmag = closed_axis_field()
    ax.plot(S, bmag, color="#00224E", lw=2)
    ratio = bmag.max() / bmag.min()
    lc = np.degrees(np.arcsin(np.sqrt(1 / ratio)))
    ax.fill_between(S, bmag, bmag.min(), alpha=0.12, color="#00224E")
    ax.set_xlabel("distance along the magnetic axis [m]")
    ax.set_ylabel("|B| [T]")
    ax.set_title(f"B · One lap of the closed loop, {S[-1]:.2f} m", loc="left",
                 fontsize=11)
    ax.annotate(f"$R_m$ = {ratio:.2f},  loss cone {lc:.1f}$^\\circ$\n"
                f"two mirror cells per lap",
                xy=(0.5, 0.42), xycoords="axes fraction", ha="center",
                fontsize=9.5, linespacing=1.4)
    ax.set_ylim(0, bmag.max() * 1.15)

    # -- C: the geometric constraint -------------------------------------
    ax = axes[1, 0]
    Ec = np.linspace(0.5, 40, 200)
    v = np.sqrt(2 * Ec * 1e3 * Q / M_D)
    rg = M_D * v / (Q * 0.1030)
    ax.plot(Ec, rg, color="#00224E", lw=2, label="gyroradius at $90^\\circ$")
    ax.plot(Ec, rg * np.sin(np.radians(lc)), color="#A69D75", lw=1.8, ls="--",
            label=f"at the loss-cone edge ({lc:.1f}$^\\circ$)")
    ax.axhline(BORE, color="#B03A2E", lw=2)
    ax.annotate(f"vessel bore {BORE:.3f} m", xy=(0.98, BORE), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", color="#B03A2E", fontsize=9)
    ax.fill_between(Ec, BORE, 0.6, color="#B03A2E", alpha=0.08)
    ax.scatter(E, [r["rg"] for r in rows], color="#00224E", zorder=5, s=28)
    ax.set_xlabel("beam energy [keV]"); ax.set_ylabel("gyroradius [m]")
    ax.set_ylim(0, 0.55)
    ax.set_title("C · The orbit has to fit the bore", loc="left", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left")

    # -- D: what that costs ----------------------------------------------
    ax = axes[1, 1]
    ax.bar(np.arange(len(E)), wall / n0, color="#B03A2E", label="vessel wall")
    ax.bar(np.arange(len(E)), (zedge + bend) / n0, bottom=wall / n0,
           color="#D9A441", label="domain edge (now zero)")
    ax.bar(np.arange(len(E)), final, bottom=(wall + zedge + bend) / n0,
           color="#00224E", label="still circulating")
    ax.set_xticks(np.arange(len(E)))
    ax.set_xticklabels([f"{e:g}" for e in E])
    ax.set_xlabel("beam energy [keV]"); ax.set_ylabel("fraction")
    ax.set_ylim(0, 1)
    ax.set_title("D · Fate after 30 $\\mu$s · the wall takes everything",
                 loc="left", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False, loc="lower left")

    for a in axes.flat:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default="racetrack.png")
    ap.add_argument("--energies", type=float, nargs="*", default=None)
    a = ap.parse_args()
    rows = analyse(a.energies if a.energies is not None else E_KEV)
    if not rows:
        raise SystemExit("no completed runs found")
    print(f"{'E[keV]':>7} {'confined':>9} {'wall':>7} {'zedge':>7} {'bend':>7} {'rg[m]':>7}")
    for r in rows:
        print(f"{r['E']:7g} {r['final']:9.3f} {r['wall']:7d} {r['zedge']:7d} "
              f"{r['bend']:7d} {r['rg']:7.3f}")
    make_figure(rows, a.plot)


if __name__ == "__main__":
    main()
