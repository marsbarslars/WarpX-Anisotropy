#!/usr/bin/env python3
"""Run and analyse the NBI injection-angle scan.

From inside the run directory:

    ../../.venv/bin/python ../../scripts/nbi_scan.py --run
    ../../.venv/bin/python ../../scripts/nbi_scan.py --plot nbi_scan.png

The sweep needs no templating — WarpX takes ParmParse overrides on the command
line, so one deck serves every angle. Each angle gets its own diagnostic prefix
so the runs do not overwrite each other's reduced diagnostics.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plasma import (confined_fraction, crossing, load_field_map,  # noqa: E402
                    load_particles, load_scraped, run_sweep, tag)

# Denser near the expected loss-cone boundary, where the interesting change is.
ANGLES = [0, 10, 20, 25, 30, 32.5, 35, 37.5, 40, 45, 50, 60, 75, 90]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--run", action="store_true", help="execute the sweep")
    p.add_argument("--plot", metavar="PATH", help="analyse and write a figure")
    p.add_argument("--inputs", default="inputs_3d_nbi_beam.txt")
    p.add_argument("--angles", type=float, nargs="*", default=None)
    p.add_argument("--max-step", type=int, default=4000)
    return p.parse_args()


def analyse(angles):
    """Per angle: confined fraction, loss channels, and the theory prediction."""
    rows = []
    fmap = None
    for a in angles:
        t = "a" + tag(a)
        d = f"diags/{t}"
        if not os.path.isdir(d):
            continue
        if fmap is None:
            from openpmd_viewer import OpenPMDTimeSeries
            fmap = load_field_map(d, OpenPMDTimeSeries(d).iterations[0])
        t_us, frac = confined_fraction(f"diags/{t}_reduced")

        scraped, where = load_scraped(f"diags/{t}_scraped", "beam", ("pitch0",))
        n_axial = n_radial = 0
        if scraped is not None and len(where):
            radial = np.isin(where, ["xlo", "xhi", "ylo", "yhi"])
            n_axial, n_radial = int((~radial).sum()), int(radial.sum())

        # Guiding-centre prediction from each particle's own birth field and
        # its pitch against the *local* B, not against the z axis.
        from openpmd_viewer import OpenPMDTimeSeries
        ts = OpenPMDTimeSeries(d)
        p0 = load_particles(d, ts.iterations[0], "beam", ("pitch0",), ts=ts)
        pred = np.nan
        if len(p0):
            b0 = fmap.interpolate(p0.x, p0.y, p0.z)
            sin2 = 1.0 - p0.pitch(b0) ** 2
            pred = float((sin2 >= np.linalg.norm(b0, axis=1)
                          / fmap.mirror_ratio()[1]).mean())

        rows.append(dict(angle=a, t_us=t_us, frac=frac,
                         final=float(frac[-1]) if frac is not None else np.nan,
                         pred=pred, axial=n_axial, radial=n_radial, path=d))
    return rows, fmap


def make_figure(rows, fmap, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theta_lc = fmap.loss_cone_deg()
    # The guiding-centre loss cone depends on which B_max a particle actually
    # reaches. On-axis is the conservative bound; a particle whose orbit wanders
    # off axis near a throat sees more field and is confined more easily.
    b_lo = fmap.mirror_ratio()[0]
    theta_off = float(np.degrees(np.arcsin(np.sqrt(b_lo / fmap.magnitude.max()))))
    ang = np.array([r["angle"] for r in rows])
    pred = np.array([r.get("pred", np.nan) for r in rows])
    final = np.array([r["final"] for r in rows])
    axial = np.array([r["axial"] for r in rows], dtype=float)
    radial = np.array([r["radial"] for r in rows], dtype=float)
    total = axial + radial
    n0 = total + final * 20000  # approximate launch count for normalisation

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    fig.suptitle(
        f"NBI angle scan · 30 keV deuterium into a mirror of ratio "
        f"{fmap.mirror_ratio()[2]:.2f}  (loss cone {theta_lc:.1f}$^\\circ$)",
        fontsize=13, y=0.98)

    # -- A: confined fraction vs angle --------------------------------------
    ax = axes[0, 0]
    ax.axvspan(theta_off, theta_lc, color="#C9A227", alpha=0.16,
               label=f"loss cone {theta_off:.1f}-{theta_lc:.1f}$^\\circ$")
    ax.axvline(theta_lc, color="#C9A227", lw=1.6)
    ax.axvline(theta_off, color="#C9A227", lw=1.6, ls=":")
    ax.plot(ang, pred, "s--", color="#A69D75", lw=1.4, ms=4,
            label="guiding centre, on-axis $B_{max}$")
    ax.plot(ang, final, "o-", color="#00224E", lw=1.8, ms=5, label="simulated")

    cross = crossing(ang, final, 0.5)
    if cross is not None:
        ax.annotate(f"50% at {cross:.1f}$^\\circ$", xy=(cross, 0.5),
                    xytext=(cross + 9, 0.42), fontsize=9, color="#00224E",
                    arrowprops=dict(arrowstyle="->", color="#00224E", lw=1))
    ax.set_xlabel("injection angle from B  [deg]")
    ax.set_ylabel("confined fraction at end of run")
    ax.set_ylim(-0.03, 1.08)
    ax.set_title("A · Confinement vs aiming angle", loc="left", fontsize=11)
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")

    # -- B: where the losses go ---------------------------------------------
    ax = axes[0, 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        f_ax = np.where(total > 0, axial / np.maximum(n0, 1), 0)
        f_rad = np.where(total > 0, radial / np.maximum(n0, 1), 0)
    ax.bar(ang, f_ax, width=2.6, color="#00224E", label="axial (mirror ends)")
    ax.bar(ang, f_rad, width=2.6, bottom=f_ax, color="#A69D75",
           label="radial (domain wall)")
    ax.axvline(theta_lc, color="#C9A227", lw=2)
    ax.set_xlabel("injection angle from B  [deg]")
    ax.set_ylabel("fraction lost")
    ax.set_title("B · Loss channel", loc="left", fontsize=11)
    ax.legend(fontsize=9, frameon=False)

    # -- C: confinement history for representative angles -------------------
    ax = axes[1, 0]
    picks = [r for r in rows
             if r["angle"] in (0, 20, 30, 35, 45, 90) and r["t_us"] is not None]
    cmap = plt.get_cmap("cividis")
    for k, r in enumerate(picks):
        ax.plot(r["t_us"], r["frac"], lw=1.8,
                color=cmap(k / max(len(picks) - 1, 1)),
                label=f"{r['angle']:g}$^\\circ$")
    ax.set_xlabel("t  [$\\mu$s]"); ax.set_ylabel("confined fraction")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title("C · Loss is prompt, not gradual", loc="left", fontsize=11)
    ax.legend(fontsize=8, frameon=False, ncol=2)

    # -- D: velocity space either side of the boundary ----------------------
    ax = axes[1, 1]
    from openpmd_viewer import OpenPMDTimeSeries
    shown = 0
    for r, colour, mark in ((_nearest(rows, 32.5), "#B03A2E", "o"),
                            (_nearest(rows, 45), "#00224E", "s")):
        if r is None:
            continue
        ts = OpenPMDTimeSeries(r["path"])
        p = load_particles(r["path"], list(ts.iterations)[-1], "beam", ts=ts)
        if not len(p):
            ax.scatter([], [], color=colour, label=f"{r['angle']:g}$^\\circ$ (all lost)")
            continue
        b = fmap.interpolate(p.x, p.y, p.z)
        vpar, vperp = p.v_par_perp(b)
        ax.scatter(vpar / 1e6, vperp / 1e6, s=2, alpha=0.25, color=colour,
                   marker=mark, linewidths=0,
                   label=f"{r['angle']:g}$^\\circ$ · {len(p)} left")
        shown += 1
    v = 3.0
    edge = np.tan(np.radians(90 - theta_lc))
    vv = np.linspace(0, v, 10)
    for sgn in (-1, 1):
        ax.plot(sgn * vv * edge, vv, c="#C9A227", lw=1.4, ls="--")
    ax.set_xlim(-v, v); ax.set_ylim(0, v)
    ax.set_xlabel("$v_\\parallel$  [10$^6$ m/s]")
    ax.set_ylabel("$v_\\perp$  [10$^6$ m/s]")
    ax.set_title("D · Survivors in velocity space", loc="left", fontsize=11)
    leg = ax.legend(fontsize=9, frameon=False, loc="upper right")
    for h in leg.legend_handles:
        try:
            h.set_alpha(1.0); h.set_sizes([24])
        except Exception:
            pass

    for a in axes.flat:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def _nearest(rows, target):
    cand = [r for r in rows if os.path.isdir(r["path"])]
    return min(cand, key=lambda r: abs(r["angle"] - target)) if cand else None


def main():
    args = parse_args()
    angles = args.angles if args.angles is not None else ANGLES

    if args.run:
        run_sweep(args.inputs, "alpha_deg", angles, species_prefix="a",
                  max_step=args.max_step, log_prefix="diags")

    if args.plot:
        rows, fmap = analyse(angles)
        if not rows:
            raise SystemExit("no completed runs found; use --run first")
        print(f"{'angle':>7} {'confined':>9} {'predicted':>10} "
              f"{'axial':>8} {'radial':>8}")
        for r in rows:
            print(f"{r['angle']:7g} {r['final']:9.3f} {r['pred']:10.3f} "
                  f"{r['axial']:8d} {r['radial']:8d}")
        b_lo, b_hi, _ = fmap.mirror_ratio()
        off = np.degrees(np.arcsin(np.sqrt(b_lo / fmap.magnitude.max())))
        print(f"\nloss cone: {fmap.loss_cone_deg():.2f} deg on-axis B_max, "
              f"{off:.2f} deg at the global B_max")
        ang = np.array([r["angle"] for r in rows])
        fin = np.array([r["final"] for r in rows])
        c = crossing(ang, fin, 0.5)
        if c is not None:
            print(f"simulated 50% crossing: {c:.2f} deg")
        make_figure(rows, fmap, args.plot)

    if not args.run and not args.plot:
        raise SystemExit("nothing to do: pass --run and/or --plot")


if __name__ == "__main__":
    main()
