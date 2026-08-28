#!/usr/bin/env python3
"""Field structure, velocity distribution and loss-cone check for a mirror run.

Run from inside the run directory:

    ../../.venv/bin/python ../../scripts/plot_mirror.py --out mirror.png

The loss-cone test is the point of the figure. A particle born where the field
is B0 escapes if its pitch angle satisfies sin^2(theta) < B0/Bmax, so the
predicted fate is known per particle from its birth position and birth pitch.
Panel F compares that prediction against what the simulation actually did.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plasma import load_field_map, load_particles  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--path", default="diags/diag1")
    p.add_argument("--species", default="ions")
    p.add_argument("--out", default="mirror.png")
    p.add_argument("--reduced", default="diags/reducedfiles")
    return p.parse_args()


def confined_series(reduced_dir):
    """(time_us, confined_fraction) from the ParticleNumber reduced diagnostic."""
    path = os.path.join(reduced_dir, "confined.txt")
    if not os.path.isfile(path):
        return None, None
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    tcol = 1 if len(header) > 1 else 0
    # the macroparticle count column, whatever it is called in this build
    cand = [i for i, h in enumerate(header)
            if "macroparticle" in h.lower() and "total" in h.lower()]
    if not cand:
        cand = [i for i, h in enumerate(header) if "macroparticle" in h.lower()]
    if not cand:
        return None, None
    n = data[:, cand[0]]
    return data[:, tcol] * 1e6, n / n[0]


def main():
    args = parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from openpmd_viewer import OpenPMDTimeSeries

    ts = OpenPMDTimeSeries(args.path)
    its = list(ts.iterations)
    first, last = its[0], its[-1]
    print(f"{len(its)} iterations, {first} -> {last}")

    field = load_field_map(args.path, first)
    b_lo, b_hi, ratio = field.mirror_ratio()
    theta_lc = field.loss_cone_deg()
    print(f"on-axis B: {b_lo:.4f} - {b_hi:.4f} T, ratio {ratio:.3f}, "
          f"loss cone {theta_lc:.2f} deg")

    p0 = load_particles(args.path, first, args.species, ("pitch0",))
    p1 = load_particles(args.path, last, args.species, ("pitch0",))
    print(f"particles: {len(p0)} -> {len(p1)}  "
          f"({100*len(p1)/len(p0):.1f}% retained)")

    b0 = field.interpolate(p0.x, p0.y, p0.z)
    b1 = field.interpolate(p1.x, p1.y, p1.z)
    vpar0, vperp0 = p0.v_par_perp(b0)
    vpar1, vperp1 = p1.v_par_perp(b1)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9))
    fig.suptitle(
        f"Magnetic mirror, {len(p0)} deuterons  |  on-axis ratio "
        f"{ratio:.2f}, loss cone {theta_lc:.1f}$^\\circ$",
        fontsize=13, y=0.98)

    # -- A: field structure -------------------------------------------------
    ax = axes[0, 0]
    z, bz = field.on_axis_bz()
    ax.plot(z, np.abs(bz), color="#00224E", lw=2)
    ax.axhline(b_lo, ls=":", c="grey", lw=1)
    ax.axhline(b_hi, ls=":", c="grey", lw=1)
    ax.fill_between(z, np.abs(bz), b_lo, alpha=0.12, color="#00224E")
    ax.set_xlabel("z [m]"); ax.set_ylabel("|B$_z$| on axis [T]")
    ax.set_title("A · Mirror field", loc="left", fontsize=11)
    ax.annotate(f"$R_m$ = {ratio:.2f}", xy=(0.5, 0.88), xycoords="axes fraction",
                ha="center", fontsize=11)

    # -- B/C: velocity space ------------------------------------------------
    vmax = np.percentile(p0.speed, 99.5) / 1e6
    bins = np.linspace(-vmax, vmax, 130)
    bins_perp = np.linspace(0, vmax, 65)

    for ax, vpar, vperp, label in (
            (axes[0, 1], vpar0, vperp0, f"B · Velocity space, step {first}"),
            (axes[0, 2], vpar1, vperp1, f"C · Velocity space, step {last}")):
        h, xe, ye = np.histogram2d(vpar / 1e6, vperp / 1e6,
                                   bins=[bins, bins_perp])
        ax.pcolormesh(xe, ye, h.T, cmap="cividis",
                      norm=LogNorm(vmin=0.5, vmax=max(h.max(), 2)))
        # loss cone: |v_par|/|v| > cos(theta_lc)
        edge = np.tan(np.radians(90 - theta_lc))
        vv = np.linspace(0, vmax, 10)
        for sgn in (-1, 1):
            ax.plot(sgn * vv * edge, vv, c="#FFEA46", lw=1.4, ls="--")
        ax.set_xlim(-vmax, vmax); ax.set_ylim(0, vmax)
        ax.set_xlabel("$v_\\parallel$ [10$^6$ m/s]")
        ax.set_ylabel("$v_\\perp$ [10$^6$ m/s]")
        ax.set_title(label, loc="left", fontsize=11)
    axes[0, 2].annotate("loss cone", xy=(0.5, 0.06), xycoords="axes fraction",
                        color="#B9A23A", ha="center", fontsize=9)

    # -- D: pitch distributions --------------------------------------------
    ax = axes[1, 0]
    a0 = np.degrees(np.arccos(np.clip(np.abs(p0.pitch(b0)), 0, 1)))
    a1 = np.degrees(np.arccos(np.clip(np.abs(p1.pitch(b1)), 0, 1)))
    pb = np.linspace(0, 90, 46)
    ax.hist(a0, bins=pb, histtype="step", lw=2, color="#575D6D",
            label=f"step {first}", density=True)
    ax.hist(a1, bins=pb, histtype="step", lw=2, color="#00224E",
            label=f"step {last}", density=True)
    ax.axvline(theta_lc, c="#C9A227", lw=2)
    ax.axvspan(0, theta_lc, color="#C9A227", alpha=0.12)
    ax.set_xlabel("pitch angle from B [deg]"); ax.set_ylabel("pdf")
    ax.set_title("D · Pitch angle: loss cone empties", loc="left", fontsize=11)
    ax.legend(fontsize=9, frameon=False)

    # -- E: confined fraction ----------------------------------------------
    ax = axes[1, 1]
    t_us, frac = confined_series(args.reduced)
    if t_us is not None:
        ax.plot(t_us, frac, color="#00224E", lw=2)
    else:
        ax.text(0.5, 0.5, "reduced diagnostic not found",
                ha="center", va="center", transform=ax.transAxes, color="grey")
    # For an isotropic population the two loss cones subtend a solid-angle
    # fraction (1 - cos(theta_lc)), so the confined fraction is cos(theta_lc).
    keep = np.cos(np.radians(theta_lc))
    ax.axhline(keep, ls="--", c="#C9A227", lw=1.6)
    ax.annotate(f"isotropic, on-axis birth: {keep:.3f}",
                xy=(0.97, keep), xycoords=("axes fraction", "data"),
                ha="right", va="bottom", fontsize=9, color="#8A6D1F")
    if frac is not None:
        ax.annotate(f"simulated: {frac[-1]:.3f}", xy=(0.97, frac[-1]),
                    xycoords=("axes fraction", "data"), ha="right", va="top",
                    fontsize=9, color="#00224E")
    ax.set_xlabel("t [$\\mu$s]"); ax.set_ylabel("confined fraction")
    ax.set_ylim(0, 1.02)
    ax.set_title("E · Confinement vs time", loc="left", fontsize=11)

    # -- F: the actual test -------------------------------------------------
    ax = axes[1, 2]
    if "pitch0" in p0.extra and "pitch0" in p1.extra:
        birth_angle = np.degrees(np.arccos(np.clip(np.abs(p0.extra["pitch0"]), 0, 1)))
        survived = np.isin(p0.ids, p1.ids)
        edges = np.linspace(0, 90, 31)
        idx = np.digitize(birth_angle, edges) - 1
        centres, rate = [], []
        for k in range(len(edges) - 1):
            m = idx == k
            if m.sum() >= 15:
                centres.append(0.5 * (edges[k] + edges[k + 1]))
                rate.append(1.0 - survived[m].mean())
        ax.plot(centres, rate, "o-", color="#00224E", ms=4, lw=1.6,
                label="simulated")

        # Per-particle criterion using each particle's own birth field rather
        # than the on-axis value: escapes if sin^2(theta) < B_birth / B_max.
        sin2 = 1.0 - np.clip(np.abs(p0.extra["pitch0"]), 0, 1) ** 2
        predicted_loss = sin2 < (np.linalg.norm(b0, axis=1) / b_hi)
        agree = (predicted_loss == ~survived).mean()
        print(f"per-particle loss-cone prediction agrees with simulation "
              f"for {100*agree:.1f}% of particles")
        print(f"  predicted lost {predicted_loss.mean():.3f}, "
              f"actually lost {1 - survived.mean():.3f}")
        prate = [predicted_loss[idx == k].mean()
                 for k in range(len(edges) - 1) if (idx == k).sum() >= 15]
        ax.plot(centres, prate, color="#A69D75", lw=1.4, ls=":",
                label="per-particle theory")
        ax.axvline(theta_lc, c="#C9A227", lw=2, label=f"theory {theta_lc:.1f}$^\\circ$")
        ax.axvspan(0, theta_lc, color="#C9A227", alpha=0.12)
        ax.set_ylim(-0.03, 1.03)
        ax.legend(fontsize=9, frameon=False, loc="center right")
    else:
        ax.text(0.5, 0.5, "pitch0 attribute missing",
                ha="center", va="center", transform=ax.transAxes, color="grey")
    ax.set_xlabel("birth pitch angle [deg]"); ax.set_ylabel("fraction lost")
    ax.set_title("F · Loss vs birth pitch", loc="left", fontsize=11)

    for a in axes.flat:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
