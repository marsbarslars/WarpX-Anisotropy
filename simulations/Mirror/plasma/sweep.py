"""Running parameter sweeps and reading their reduced diagnostics.

WarpX takes ParmParse overrides on the command line, so a sweep needs no
templating — one deck serves every point. The only discipline required is giving
each point its own output prefixes: `diag1.file_prefix`, `scraped.file_prefix`
*and* `reduced_diags.path`. Forgetting the last silently overwrites the previous
point's `confined.txt`, since reduced diagnostics default to `diags/reducedfiles/`
regardless of where the other diagnostics went.
"""

from __future__ import annotations

import os
import subprocess

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
WARPX = os.path.join(_HERE, os.pardir, "warpx")


def tag(value) -> str:
    """Filesystem-safe label for a sweep point ('32.5' -> '32p5')."""
    return f"{value:g}".replace(".", "p").replace("-", "m")


def run_sweep(inputs, param, values, species_prefix="", extra=(), max_step=None,
              log_prefix="sweep"):
    """Run `inputs` once per value of `param`, each into its own diagnostics."""
    for i, v in enumerate(values, 1):
        t = f"{species_prefix}{tag(v)}"
        print(f"[{i}/{len(values)}] {param} = {v} -> diags/{t}", flush=True)
        cmd = [WARPX, "run", inputs, f"{param}={v}",
               f"diag1.file_prefix=diags/{t}",
               f"scraped.file_prefix=diags/{t}_scraped",
               f"reduced_diags.path=./diags/{t}_reduced/"]
        if max_step is not None:
            cmd.append(f"max_step={max_step}")
        cmd += list(extra)
        with open(f"{log_prefix}_{t}.log", "w") as log:
            r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            raise SystemExit(f"{param}={v} failed; see {log_prefix}_{t}.log")


def confined_fraction(reduced_dir):
    """(time_us, fraction_remaining) from a ParticleNumber reduced diagnostic."""
    path = os.path.join(reduced_dir, "confined.txt")
    if not os.path.isfile(path):
        return None, None
    with open(path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    col = next((i for i, h in enumerate(header)
                if "macroparticle" in h.lower() and "total" in h.lower()), None)
    if col is None:
        return None, None
    n = data[:, col]
    return data[:, 1] * 1e6, n / n[0]


def crossing(x, y, level):
    """Linear interpolation of where y first rises through `level`."""
    x, y = np.asarray(x), np.asarray(y)
    for i in range(len(x) - 1):
        if y[i] < level <= y[i + 1]:
            f = (level - y[i]) / (y[i + 1] - y[i])
            return float(x[i] + f * (x[i + 1] - x[i]))
    return None
