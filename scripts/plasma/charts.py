"""Live chart panel composited beside the 3D render.

Two panels:

* a **static** |B| scale bar — the legend for how the field lines in the 3D view
  are coloured, so the two read as one picture;
* the (v_par, v_perp) distribution over *every* particle in the frame, where the
  loss cone empties out over the run.

Colour is doing two different jobs here, so it uses two maps rather than
overloading one. Field strength stays on cividis, matching the field lines it
labels. Particle density uses an inferno-like ramp that starts at the panel
background, so empty velocity space reads as empty rather than as a low value.
The loss cone is drawn in a cool accent that neither map contains, so it never
looks like part of the data.

Rendered with Agg into an RGB array so it can be stacked next to a PyVista
screenshot; nothing here opens a window.
"""

from __future__ import annotations

import numpy as np

# Panel palette. GROUND matches the 3D view's background so the composite reads
# as one canvas; the rest is a cool grey ramp biased toward it.
GROUND = "#07090F"
PANEL = "#0D1119"
RULE = "#2A3446"
INK = "#E6EAF2"
INK_DIM = "#8B96AA"
ACCENT = "#5FD3F3"      # loss cone: absent from both colormaps


class ChartPanel:
    """Fixed-size matplotlib panel, redrawn per frame.

    Axis limits and the density colour ceiling are locked at construction:
    left to autoscale they drift as particles are lost, and the charts jitter.
    """

    def __init__(self, width_px, height_px, b_range, v_scale, loss_cone_deg,
                 dpi=100):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        self.dpi = dpi
        self.b_lo, self.b_hi = b_range
        self.v_scale = v_scale
        self.theta_lc = loss_cone_deg
        self._ceiling = None

        self.fig = plt.figure(
            figsize=(width_px / dpi, height_px / dpi), dpi=dpi,
            facecolor=GROUND)
        self.fig.subplots_adjust(left=0.17, right=0.94, top=0.92, bottom=0.07,
                                 hspace=0.42)

        gs = self.fig.add_gridspec(2, 1, height_ratios=[0.5, 3.9])
        self.ax_b = self.fig.add_subplot(gs[0])
        self.ax_v = self.fig.add_subplot(gs[1])

        for ax in (self.ax_b, self.ax_v):
            ax.set_facecolor(PANEL)
            for spine in ax.spines.values():
                spine.set_color(RULE)
            ax.tick_params(colors=INK_DIM, labelsize=8, length=3)
            ax.xaxis.label.set_color(INK_DIM)
            ax.yaxis.label.set_color(INK_DIM)

        self._setup_b_bar()
        self._setup_velocity()

    # ------------------------------------------------------ |B| bar (static)

    def _setup_b_bar(self):
        ax = self.ax_b
        grad = np.linspace(0, 1, 256)[None, :]
        ax.imshow(grad, aspect="auto", cmap="cividis",
                  extent=[self.b_lo, self.b_hi, 0, 1])
        ax.set_yticks([])
        ax.set_xlim(self.b_lo, self.b_hi)
        ax.set_title("Magnetic field strength", fontsize=11, loc="left",
                     pad=24, fontweight="medium", color=INK)
        _subtitle(ax, "field-line colouring at left · |B| in tesla")

        # Endpoints go *below* the bar: white-on-yellow is illegible at the
        # throat end, and the two values are what set the loss cone.
        for value, label, ha in ((self.b_lo, "midplane", "left"),
                                 (self.b_hi, "throat", "right")):
            ax.annotate(label,
                        xy=(value, 0), xycoords=("data", "axes fraction"),
                        xytext=(0, -21), textcoords="offset points",
                        ha=ha, va="top", fontsize=8.5, color=INK)

    # ------------------------------------------------------- velocity space

    def _setup_velocity(self):
        import matplotlib.colors as mcolors

        ax = self.ax_v
        v = self.v_scale
        # Coarser than looks natural on purpose. At 110x56 there are ~4 counts
        # per occupied bin, so neighbouring bins differ by Poisson noise of the
        # same size as the signal, and a linear colour ramp turns each integer
        # count into its own visible band.
        self.xedges = np.linspace(-v, v, 78)
        self.yedges = np.linspace(0, v, 40)

        # Start the ramp exactly at the panel colour so empty bins disappear.
        base = __import__("matplotlib").colormaps["inferno"](np.linspace(0, 1, 256))
        base[0] = mcolors.to_rgba(PANEL)
        self._cmap = mcolors.ListedColormap(base)

        self._v_img = ax.pcolormesh(
            self.xedges, self.yedges,
            np.zeros((len(self.yedges) - 1, len(self.xedges) - 1)),
            cmap=self._cmap, vmin=0, vmax=1)

        edge = np.tan(np.radians(90 - self.theta_lc))
        vv = np.linspace(0, v, 10)
        for sgn in (-1, 1):
            ax.plot(sgn * vv * edge, vv, c=ACCENT, lw=1.3, ls=(0, (5, 4)),
                    alpha=0.95, zorder=5)

        ax.set_xlim(-v, v); ax.set_ylim(0, v)
        ax.set_xlabel("$v_\\parallel$  [10$^6$ m/s]", fontsize=9, labelpad=2)
        ax.set_ylabel("$v_\\perp$  [10$^6$ m/s]", fontsize=9, labelpad=2)
        ax.set_title("Velocity distribution", fontsize=11, loc="left", pad=24,
                     fontweight="medium", color=INK)
        _subtitle(ax, f"all particles · loss cone {self.theta_lc:.1f}$^\\circ$",
                  color=ACCENT)
        # The histogram's bright end can sit under this text, so give it a
        # ground of its own rather than relying on what happens to be behind.
        self._v_text = ax.text(0.97, 0.95, "", transform=ax.transAxes,
                               ha="right", va="top", color=INK, fontsize=9,
                               zorder=6,
                               bbox=dict(facecolor=GROUND, edgecolor="none",
                                         alpha=0.75, pad=3.5))

    def _density(self, vpar, vperp):
        """Smoothed 2D density. Raw counts are a noisy estimator at this
        sample size; a light Gaussian makes the panel show the distribution
        rather than the shot noise on it."""
        h, _, _ = np.histogram2d(vpar / 1e6, vperp / 1e6,
                                 bins=[self.xedges, self.yedges])
        try:
            from scipy.ndimage import gaussian_filter
        except ImportError:
            return h
        return gaussian_filter(h, sigma=1.1, mode="nearest")

    def calibrate(self, samples):
        """Fix the colour ceiling from frames spread across the run.

        Locking it on frame 0 clips badly: the population starts spread out and
        then concentrates, so the final peak is several times the initial one
        and most of the distribution saturates to flat white.
        """
        tops = []
        for vpar, vperp in samples:
            if len(vpar) == 0:
                continue
            d = self._density(vpar, vperp)
            if np.any(d > 0):
                tops.append(float(np.percentile(d[d > 0], 99)))
        self._ceiling = max(tops) if tops else 1.0
        self._v_img.set_clim(0, self._ceiling)

    def _update_velocity(self, vpar, vperp, n_total):
        d = self._density(vpar, vperp)
        if self._ceiling is None:      # not calibrated; fall back to this frame
            self._ceiling = max(float(np.percentile(d[d > 0], 99))
                                if np.any(d > 0) else 1.0, 1.0)
            self._v_img.set_clim(0, self._ceiling)
        self._v_img.set_array(d.T.ravel())
        self._v_text.set_text(f"{len(vpar):,} of {n_total:,} confined")

    # ------------------------------------------------------------------ draw

    def render(self, bmag, vpar, vperp, n_total) -> np.ndarray:
        """Update the velocity panel and return the whole panel as RGB.

        ``bmag`` is accepted and ignored — the field bar is a static legend.
        """
        self._update_velocity(vpar, vperp, n_total)
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba())[:, :, :3].copy()

    def close(self):
        import matplotlib.pyplot as plt
        plt.close(self.fig)


def _subtitle(ax, text, color=INK_DIM):
    """Caption under a title, offset in points so panel height cannot move it."""
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction",
                xytext=(0, 5), textcoords="offset points",
                ha="left", va="bottom", color=color, fontsize=9)


def stack_side_by_side(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Join two RGB frames horizontally, padding to the taller one."""
    h = max(left.shape[0], right.shape[0])
    out = np.zeros((h, left.shape[1] + right.shape[1], 3), dtype=np.uint8)
    out[:left.shape[0], :left.shape[1]] = left
    out[:right.shape[0], left.shape[1]:] = right
    return out
