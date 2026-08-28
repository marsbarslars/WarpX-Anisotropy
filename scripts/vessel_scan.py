#!/usr/bin/env python3
"""Vessel-radius scan: what a real wall costs, and where the load lands.

From inside the run directory:

    ../../.venv/bin/python ../../scripts/vessel_scan.py --run
    ../../.venv/bin/python ../../scripts/vessel_scan.py --plot vessel.png

Phase 1 of docs/nbi-plan.md. Until now the domain edge stood in for a wall, so
radial losses were an artefact of the box. With an embedded boundary the losses
split cleanly: axial ones leave through the mirror throats, radial ones land on
the vessel and carry a surface normal, which is what makes a wall-load map
possible.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plasma import (confined_fraction, load_field_map, load_particles,  # noqa: E402
                    load_scraped, run_sweep, tag)

RADII = [0.25, 0.35, 0.45, 0.55, 0.70, 0.90]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run", action="store_true")
    p.add_argument("--plot", metavar="PATH")
    p.add_argument("--inputs", default="inputs_3d_vessel.txt")
    p.add_argument("--radii", type=float, nargs="*", default=None)
    p.add_argument("--max-step", type=int, default=4000)
    return p.parse_args()


def analyse(radii):
    rows, fmap = [], None
    for r in radii:
        t = tag(r)
        d = f"diags/{t}"
        if not os.path.isdir(d):
            continue
        if fmap is None:
            from openpmd_viewer import OpenPMDTimeSeries
            fmap = load_field_map(d, OpenPMDTimeSeries(d).iterations[0])
        t_us, frac = confined_fraction(f"diags/{t}_reduced")
        # Particles born outside the vessel are deleted at initialisation and
        # never reach the scraping buffer, so the launched count is whatever
        # survived setup - not the requested npart.
        n_launched = _launched(f"diags/{t}_reduced")
        scraped, where = load_scraped(f"diags/{t}_scraped", "ions",
                                      ("pitch0", "nx", "ny", "nz"))
        n_eb = n_axial = n_box = 0
        eb = None
        if scraped is not None and len(where):
            m_eb = where == "eb"
            m_ax = np.isin(where, ["zlo", "zhi"])
            n_eb, n_axial = int(m_eb.sum()), int(m_ax.sum())
            n_box = int(len(where) - n_eb - n_axial)
            if m_eb.any():
                eb = dict(x=scraped.x[m_eb], y=scraped.y[m_eb],
                          z=scraped.z[m_eb],
                          pitch0=scraped.extra.get("pitch0",
                                                   np.zeros(len(scraped)))[m_eb])
        rows.append(dict(radius=r, t_us=t_us, frac=frac, n0=n_launched,
                         final=float(frac[-1]) if frac is not None else np.nan,
                         eb=n_eb, axial=n_axial, box=n_box, eb_hits=eb, path=d))
    return rows, fmap


def _launched(reduced_dir):
    """Macroparticles actually present at step 0, after any born inside the
    vessel wall have been culled."""
    path = os.path.join(reduced_dir, "confined.txt")
    if not os.path.isfile(path):
        return np.nan
    d = np.loadtxt(path, skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    return int(d[0, 2])


def make_figure(rows, fmap, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta_lc = fmap.loss_cone_deg()
    rad = np.array([r["radius"] for r in rows])
    final = np.array([r["final"] for r in rows])
    eb = np.array([r["eb"] for r in rows], float)
    axial = np.array([r["axial"] for r in rows], float)
    n0 = np.array([r["n0"] for r in rows], float)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"Vacuum vessel scan · isotropic 30 keV deuterons, mirror ratio "
        f"{fmap.mirror_ratio()[2]:.2f} (loss cone {theta_lc:.1f}$^\\circ$)",
        fontsize=13, y=0.98)

    # -- A: confinement vs vessel radius ------------------------------------
    ax = axes[0, 0]
    ax.plot(rad, final, "o-", color="#00224E", lw=1.8, ms=5, label="confined")
    keep = np.cos(np.radians(theta_lc))
    ax.axhline(keep, ls="--", c="#C9A227", lw=1.6)
    ax.annotate(f"mirror-only limit {keep:.3f}", xy=(rad[-1], keep),
                xytext=(0, 5), textcoords="offset points", ha="right",
                fontsize=9, color="#8A6D1F")
    ax.set_xlabel("vessel radius  [m]"); ax.set_ylabel("confined fraction")
    ax.set_ylim(0, 1.02)
    ax.set_title("A · The wall costs confinement", loc="left", fontsize=11)
    ax.legend(fontsize=9, frameon=False, loc="lower right")

    # -- B: loss channel ----------------------------------------------------
    ax = axes[0, 1]
    w = 0.045
    ax.bar(rad, axial / n0, width=w, color="#00224E", label="axial (throats)")
    ax.bar(rad, eb / n0, width=w, bottom=axial / n0, color="#B03A2E",
           label="vessel wall")
    ax.set_xlabel("vessel radius  [m]"); ax.set_ylabel("fraction lost")
    ax.set_title("B · Loss channel", loc="left", fontsize=11)
    ax.legend(fontsize=9, frameon=False)

    # -- C: wall load map ---------------------------------------------------
    ax = axes[1, 0]
    mid = min(rows, key=lambda r: abs(r["radius"] - 0.45))
    _ = mid
    hits = mid["eb_hits"]
    if hits is not None and len(hits["z"]):
        phi = np.degrees(np.arctan2(hits["y"], hits["x"]))
        hh, xe, ye = np.histogram2d(hits["z"], phi, bins=[48, 36])
        pcm = ax.pcolormesh(xe, ye, hh.T, cmap="inferno")
        fig.colorbar(pcm, ax=ax, label="hits per cell", pad=0.02)
        ax.set_ylabel("azimuth  [deg]")
    ax.set_xlabel("z  [m]")
    ax.set_title(f"C · Wall load, r = {mid['radius']:g} m", loc="left",
                 fontsize=11)

    # -- D: axial profile of the wall load ----------------------------------
    ax = axes[1, 1]
    zg, bzg = fmap.on_axis_bz()
    for r in rows:
        h = r["eb_hits"]
        if h is None or not len(h["z"]):
            continue
        cnt, edges = np.histogram(h["z"], bins=40, range=(zg[0], zg[-1]))
        ax.plot(0.5 * (edges[1:] + edges[:-1]), cnt / max(cnt.sum(), 1),
                lw=1.6, label=f"r = {r['radius']:g} m")
    ax.set_xlabel("z  [m]"); ax.set_ylabel("fraction of wall hits")
    ax.set_title("D · Where the wall load lands", loc="left", fontsize=11)
    ax.legend(fontsize=8, frameon=False, ncol=2)
    axb = ax.twinx()
    axb.plot(zg, np.abs(bzg), color="#A69D75", lw=1.2, ls=":")
    axb.set_ylabel("|B$_z$| on axis  [T]", color="#8A7B4A", fontsize=9)
    axb.tick_params(axis="y", colors="#8A7B4A", labelsize=8)
    axb.spines["top"].set_visible(False)

    for a in axes.flat:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def main():
    args = parse_args()
    radii = args.radii if args.radii is not None else RADII

    if args.run:
        # rb defaults to ra in the deck, so overriding ra alone keeps the
        # vessel circular across the sweep.
        run_sweep(args.inputs, "ra", radii, max_step=args.max_step,
                  log_prefix="vessel")

    if args.plot:
        rows, fmap = analyse(radii)
        if not rows:
            raise SystemExit("no completed runs found; use --run first")
        print(f"{'radius':>7} {'launched':>9} {'confined':>9} {'axial':>7} "
              f"{'wall':>7} {'box':>5} {'balance':>8}")
        for r in rows:
            n_conf = int(round(r["final"] * r["n0"]))
            bal = n_conf + r["axial"] + r["eb"] + r["box"] - r["n0"]
            print(f"{r['radius']:7g} {r['n0']:9d} {r['final']:9.3f} "
                  f"{r['axial']:7d} {r['eb']:7d} {r['box']:5d} {bal:+8d}")
        make_figure(rows, fmap, args.plot)

    if not args.run and not args.plot:
        raise SystemExit("nothing to do: pass --run and/or --plot")


if __name__ == "__main__":
    main()
