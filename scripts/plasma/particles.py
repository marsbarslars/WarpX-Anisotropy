"""Loading WarpX particle output and deriving velocity-space quantities.

Two units traps live here, both of which produce plausible-looking but wrong
numbers if ignored:

* openPMD stores ``momentum`` as relativistic momentum in kg m/s, so velocity is
  ``p / (gamma m)``. At NBI energies gamma is within 1e-4 of 1, but the division
  is cheap and keeps the code honest at higher energy.
* The parser that fills custom attributes at injection reports ``ux, uy, uz`` in
  WarpX's internal units of gamma*v (m/s), *not* the gamma*beta of the
  documentation. So a ``speed0`` attribute comes back in m/s.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

C = 2.99792458e8


@dataclass
class Particles:
    """Particle state at one iteration, in SI units."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray
    ids: np.ndarray
    mass: float
    extra: dict

    def __len__(self) -> int:
        return len(self.x)

    @property
    def speed(self) -> np.ndarray:
        return np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

    @property
    def energy_kev(self) -> np.ndarray:
        return 0.5 * self.mass * self.speed**2 / 1.602176634e-19 / 1e3

    def pitch(self, bfield: np.ndarray) -> np.ndarray:
        """cos(pitch angle) against a per-particle B, shape (N, 3).

        +1 is field-aligned, 0 is purely perpendicular. This is the quantity the
        loss cone is defined on.
        """
        bmag = np.linalg.norm(bfield, axis=1)
        vpar = (self.vx * bfield[:, 0]
                + self.vy * bfield[:, 1]
                + self.vz * bfield[:, 2]) / np.where(bmag > 0, bmag, np.inf)
        return vpar / np.where(self.speed > 0, self.speed, np.inf)

    def v_par_perp(self, bfield: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(v_parallel, v_perpendicular) — the axes of the canonical fast-ion plot."""
        bmag = np.linalg.norm(bfield, axis=1)
        safe = np.where(bmag > 0, bmag, np.inf)
        vpar = (self.vx * bfield[:, 0]
                + self.vy * bfield[:, 1]
                + self.vz * bfield[:, 2]) / safe
        vperp = np.sqrt(np.maximum(self.speed**2 - vpar**2, 0.0))
        return vpar, vperp


def _velocity_from_momentum(series, species, iteration, mass):
    px, py, pz = (series.get_particle([f"u{c}"], species=species,
                                      iteration=iteration)[0] for c in "xyz")
    # openpmd-viewer returns momentum/(m c) for ux; convert to velocity.
    gamma = np.sqrt(1.0 + px**2 + py**2 + pz**2)
    return (px * C / gamma, py * C / gamma, pz * C / gamma)


def load_particles(series_path: str, iteration: int, species: str = "ions",
                   extra_records=(), ts=None) -> Particles:
    """Load one iteration of a Full diagnostic.

    Pass an existing ``ts`` when looping over many iterations — constructing an
    OpenPMDTimeSeries rescans the directory, which dominates the cost of a
    per-frame read.
    """
    from openpmd_viewer import OpenPMDTimeSeries

    if ts is None:
        ts = OpenPMDTimeSeries(series_path)
    x, y, z, ids = ts.get_particle(["x", "y", "z", "id"],
                                   species=species, iteration=iteration)

    # A species can be entirely absorbed — every particle inside the loss cone,
    # for instance — leaving a valid but empty record. Return an empty set
    # rather than indexing into nothing.
    if len(x) == 0:
        e = np.array([])
        return Particles(e, e, e, e, e, e, e.astype(int), float("nan"),
                         {r: e for r in extra_records})

    mass = ts.get_particle(["mass"], species=species, iteration=iteration)[0][0]
    vx, vy, vz = _velocity_from_momentum(ts, species, iteration, mass)

    extra = {}
    for rec in extra_records:
        try:
            extra[rec] = ts.get_particle([rec], species=species,
                                         iteration=iteration)[0]
        except Exception:
            pass  # attribute absent in this diagnostic

    return Particles(x, y, z, vx, vy, vz, ids, float(mass), extra)


def load_scraped(diag_dir: str, species: str = "ions", extra_records=()):
    """Load every particle a BoundaryScraping diagnostic caught, merged across
    boundaries. Returns (Particles, boundary_label array)."""
    import os
    from openpmd_viewer import OpenPMDTimeSeries

    chunks, labels = [], []
    for name in sorted(os.listdir(diag_dir)):
        sub = os.path.join(diag_dir, name)
        if not name.startswith("particles_at_") or not os.path.isdir(sub):
            continue
        try:
            ts = OpenPMDTimeSeries(sub)
        except Exception:
            continue
        for it in ts.iterations:
            try:
                p = load_particles(sub, it, species, extra_records)
            except Exception:
                continue
            if len(p):
                chunks.append(p)
                labels.append(np.full(len(p), name[len("particles_at_"):]))

    if not chunks:
        return None, np.array([])

    cat = lambda attr: np.concatenate([getattr(c, attr) for c in chunks])
    extra = {k: np.concatenate([c.extra[k] for c in chunks])
             for k in chunks[0].extra}
    merged = Particles(
        cat("x"), cat("y"), cat("z"),
        cat("vx"), cat("vy"), cat("vz"), cat("ids"),
        chunks[0].mass, extra,
    )
    return merged, np.concatenate(labels)
