#!/usr/bin/env python3
"""Plot the matched mirror campaign and a near-loss-cone velocity-space case."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


HERE = Path(__file__).resolve().parent
OUT = HERE / "visualizations"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: str) -> float:
    return float(value) if value not in {"", "None"} else np.nan


def plot_summary() -> None:
    rows = read_csv(HERE / "campaign_comparison.csv")
    angles = np.array([number(row["angle_deg"]) for row in rows])
    loss_cone = 34.548
    control_survival = np.array([number(row["collisionless_survival"]) for row in rows])
    collision_survival = np.array([number(row["collisional_survival"]) for row in rows])
    control_energy = np.array([number(row["collisionless_confined_energy_fraction"]) for row in rows])
    collision_energy = np.array([number(row["collisional_confined_energy_fraction"]) for row in rows])
    control_pitch = np.array([number(row["collisionless_pitch_variance"]) for row in rows])
    collision_pitch = np.array([number(row["collisional_pitch_variance"]) for row in rows])

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    style = {"marker": "o", "linewidth": 2.2}
    for ax in axes.flat:
        ax.axvline(loss_cone, color="#777777", linestyle=":", linewidth=1.6)
        ax.grid(alpha=0.22)
        ax.set_xlabel("Injection angle relative to local B [deg]")

    axes[0, 0].plot(angles, 100 * control_survival, color="#2878b5", label="collisionless", **style)
    axes[0, 0].plot(angles, 100 * collision_survival, color="#e26a2c", label="collisional", **style)
    axes[0, 0].set_ylabel("Beam survival at 299.2 µs [%]")
    axes[0, 0].set_ylim(-3, 104)
    axes[0, 0].legend(frameon=False)

    axes[0, 1].plot(angles, 100 * control_energy, color="#2878b5", **style)
    axes[0, 1].plot(angles, 100 * collision_energy, color="#e26a2c", **style)
    axes[0, 1].set_ylabel("Confined beam energy / initial energy [%]")
    axes[0, 1].set_ylim(-3, 104)

    axes[1, 0].plot(angles, control_pitch, color="#2878b5", **style)
    axes[1, 0].plot(angles, collision_pitch, color="#e26a2c", **style)
    axes[1, 0].set_ylabel("Survivor pitch variance, Var($v_\\parallel/v$)")

    delta = 100 * (collision_survival - control_survival)
    axes[1, 1].axhline(0.0, color="black", linewidth=1.0)
    axes[1, 1].bar(angles, delta, width=6.5, color=np.where(delta >= 0, "#3a9d5d", "#bb4d58"))
    axes[1, 1].set_ylabel("Collision-associated survival change [percentage points]")

    fig.suptitle(
        "Magnetic-mirror injection-angle campaign · 5,000 ions/case · $n_i=10^{16}$ m$^{-3}$\n"
        "Dotted line: theoretical 34.55° loss-cone boundary",
        fontsize=14,
    )
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "campaign_summary.png", dpi=180)
    plt.close(fig)


def plot_presentation_summary() -> None:
    """Render only the two observables needed for the two-minute presentation."""
    rows = read_csv(HERE / "campaign_comparison.csv")
    angles = np.array([number(row["angle_deg"]) for row in rows])
    loss_cone = 34.548
    control_survival = np.array([number(row["collisionless_survival"]) for row in rows])
    collision_survival = np.array([number(row["collisional_survival"]) for row in rows])
    control_energy = np.array([number(row["collisionless_confined_energy_fraction"]) for row in rows])
    collision_energy = np.array([number(row["collisional_confined_energy_fraction"]) for row in rows])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    style = {"marker": "o", "linewidth": 2.6, "markersize": 7}
    for ax in axes:
        ax.axvline(loss_cone, color="#666666", linestyle=":", linewidth=2.0)
        ax.text(
            loss_cone + 0.9,
            4.0,
            "34.55° loss cone",
            rotation=90,
            color="#555555",
            va="bottom",
            fontsize=9,
        )
        ax.grid(alpha=0.22)
        ax.set_xlabel("Injection angle relative to local B [deg]")
        ax.set_ylim(-3, 104)

    axes[0].plot(angles, 100 * control_survival, color="#2878b5", label="collisionless", **style)
    axes[0].plot(angles, 100 * collision_survival, color="#e26a2c", label="collisional", **style)
    axes[0].set_ylabel("Beam survival at 299.2 µs [%]")
    axes[0].set_title("How injection angle controls confinement")
    axes[0].legend(frameon=False, loc="lower right")
    axes[0].annotate(
        "collisions rescue\nsome passing ions",
        xy=(30, 33.46),
        xytext=(17, 48),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )
    axes[0].annotate(
        "collisions scatter some\ntrapped ions into loss",
        xy=(40, 66.22),
        xytext=(49, 52),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )

    axes[1].plot(angles, 100 * control_energy, color="#2878b5", label="collisionless", **style)
    axes[1].plot(angles, 100 * collision_energy, color="#e26a2c", label="collisional", **style)
    axes[1].set_ylabel("Beam energy remaining in domain / initial [%]")
    axes[1].set_title("Loss plus collisional slowing")
    axes[1].legend(frameon=False, loc="lower right")
    axes[1].annotate(
        "collisionless energy follows survival:\nescaping ions carry their energy out",
        xy=(50, 98.16),
        xytext=(39, 76),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )
    axes[1].annotate(
        "collisional beam also slows down",
        xy=(60, 57.59),
        xytext=(45, 35),
        arrowprops={"arrowstyle": "->", "color": "#555555"},
        fontsize=9,
    )

    fig.suptitle(
        "Magnetic-mirror angle scan · 5,000 ions/case · $n_i=10^{16}$ m$^{-3}$",
        fontsize=15,
    )
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "campaign_summary_presentation.png", dpi=180)
    plt.close(fig)


def plot_velocity_space() -> None:
    control_case = HERE / "runs" / "campaign_angle_035_collisionless"
    collision_case = HERE / "runs" / "campaign_angle_035_collisional"
    control = np.load(control_case / "velocity_space.npz")
    collision = np.load(collision_case / "velocity_space.npz")
    control_summary = json.loads((control_case / "campaign_summary.json").read_text(encoding="utf-8"))
    collision_summary = json.loads((collision_case / "campaign_summary.json").read_text(encoding="utf-8"))
    mirror_ratio = float(control_summary["field_mirror_ratio"])

    panels = [
        (control["vparallel_vperp_counts"][0], "Initial beam", 5000),
        (
            control["vparallel_vperp_counts"][-1],
            "Collisionless · 299.2 µs",
            control_summary["final_beam_count"],
        ),
        (
            collision["vparallel_vperp_counts"][-1],
            "Collisional · 299.2 µs",
            collision_summary["final_beam_count"],
        ),
    ]
    x_edges = control["vparallel_edges_m_s"] / 1000.0
    y_edges = control["vperp_edges_m_s"] / 1000.0
    probabilities = [counts.T / max(np.sum(counts), 1) for counts, _, _ in panels]
    positive = np.concatenate([values[values > 0] for values in probabilities])
    norm = LogNorm(vmin=max(float(positive.min()), 1.0e-5), vmax=float(max(v.max() for v in probabilities)))
    loss_boundary_vparallel = np.linspace(0.0, float(np.max(np.abs(x_edges))), 240)
    loss_boundary_vperp = loss_boundary_vparallel / np.sqrt(mirror_ratio - 1.0)

    fig = plt.figure(figsize=(12, 7), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=(3, 1.35))
    images = []
    for column, ((_, title, count), probability) in enumerate(zip(panels, probabilities)):
        ax = fig.add_subplot(grid[0, column])
        images.append(ax.pcolormesh(x_edges, y_edges, probability, shading="auto", cmap="magma", norm=norm))
        ax.set_title(f"{title}\nN = {count}")
        ax.set_xlabel("$v_\\parallel$ [km/s]")
        if column == 0:
            ax.set_ylabel("$v_\\perp$ [km/s]")
        for sign in (-1.0, 1.0):
            ax.plot(
                sign * loss_boundary_vparallel,
                loss_boundary_vperp,
                color="#777777",
                linestyle=":",
                linewidth=1.8,
                label="loss-cone boundary" if column == 0 and sign > 0 else None,
            )
        if column == 0:
            ax.legend(frameon=False, loc="upper left")
        ax.grid(alpha=0.12)
    fig.colorbar(images[-1], ax=fig.axes[:3], label="Probability per velocity bin", shrink=0.88)

    ax = fig.add_subplot(grid[1, :])
    xi_edges = control["xi_edges"]
    xi_center = 0.5 * (xi_edges[:-1] + xi_edges[1:])
    for archive, index, label, color in (
        (control, 0, "initial", "#555555"),
        (control, -1, "collisionless final", "#2878b5"),
        (collision, -1, "collisional final", "#e26a2c"),
    ):
        counts = archive["xi_counts"][index].astype(float)
        ax.step(xi_center, counts / max(counts.sum(), 1.0), where="mid", label=label, color=color, linewidth=2)
    xi_boundary = np.cos(np.deg2rad(34.548))
    ax.axvline(xi_boundary, color="#777777", linestyle=":")
    ax.axvline(-xi_boundary, color="#777777", linestyle=":")
    ax.set_xlabel("Local pitch cosine $\\xi=v_\\parallel/v$")
    ax.set_ylabel("Probability")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, ncol=3)
    fig.suptitle("35° near-loss-cone velocity-space redistribution")
    fig.savefig(OUT / "campaign_velocity_space_035.png", dpi=180)
    plt.close(fig)


def main() -> None:
    plot_summary()
    plot_presentation_summary()
    plot_velocity_space()
    print(OUT / "campaign_summary.png")
    print(OUT / "campaign_summary_presentation.png")
    print(OUT / "campaign_velocity_space_035.png")


if __name__ == "__main__":
    main()
