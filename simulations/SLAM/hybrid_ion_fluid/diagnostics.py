"""Cell-defined ion diagnostics. Parallel/perpendicular refer to local total B.

Temperatures use physical weights and CELL-local mean-flow subtraction.
This empirical pressure tensor is not a Maxwellian fit, nor an unbiased
finite-sample estimator. Low-occupancy cells are explicitly reported.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

QE, MP, C, KB = 1.602176634e-19, 1.67262192595e-27, 299792458.0, 1.380649e-23


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


@dataclass(frozen=True)
class Grid:
    lo: tuple = (-1., -1., 0.)
    hi: tuple = (1., 1., 5.)
    shape: tuple = (32, 32, 80)

    @property
    def dx(self):
        return (np.array(self.hi)-self.lo)/self.shape

    @property
    def volume(self):
        return float(np.prod(self.dx))

    def centers(self):
        ijk = np.indices(self.shape).reshape(3, -1).T
        return np.array(self.lo)+(ijk+.5)*self.dx

    def ids(self, xyz):
        ijk = np.floor((xyz-np.array(self.lo))/self.dx).astype(int)
        valid = np.all((ijk >= 0) & (ijk < np.array(self.shape)), axis=1)
        ids = np.full(len(xyz), -1, dtype=int)
        ids[valid] = np.ravel_multi_index(ijk[valid].T, self.shape)
        return ids

    def as_dict(self):
        return dict(lo_m=list(self.lo), hi_m=list(self.hi), shape=list(self.shape),
                    dx_m=self.dx.tolist(), indexing="zero-based i=x, j=y, k=z; upper bounds exclusive")


def regions(grid):
    if grid != Grid():
        raise ValueError("ROI indices must be redesigned for a different mesh/vessel")
    ijk = np.indices(grid.shape).reshape(3, -1).T
    x, y, z = grid.centers().T
    center = (x*x+y*y < .25**2) & (ijk[:, 2] >= 36) & (ijk[:, 2] < 44)
    injection = np.all((ijk >= (12, 14, 14)) & (ijk < (16, 18, 18)), axis=1)
    masks = dict(center=center, near_injection=injection)
    reasons = dict(center="Central low-B bulk, separated from source; compare with old central probe",
                   near_injection="Contains x=-0.2 source and its y/z aperture, includes immediate downstream cells")
    result = {}
    for name, mask in masks.items():
        cells = ijk[mask]
        result[name] = dict(cell_ids=np.flatnonzero(mask).tolist(), ijk=cells.tolist(),
                            count=len(cells), volume_m3=len(cells)*grid.volume,
                            bounds_m=[(np.array(grid.lo)+cells.min(0)*grid.dx).tolist(),
                                      (np.array(grid.lo)+(cells.max(0)+1)*grid.dx).tolist()],
                            reason=reasons[name], selection="WHOLE Cartesian cells; no fractional cylinder volume")
    return result


def moments(xyz, velocity, weights, Bcells, grid, region):
    xyz, velocity, weights = map(lambda x: np.asarray(x, float), (xyz, velocity, weights))
    n = int(np.prod(grid.shape))
    if xyz.shape != velocity.shape or xyz.shape != (len(weights), 3) or Bcells.shape != (n, 3):
        raise ValueError("Particle or B-field array shape mismatch")
    if not all(np.isfinite(a).all() for a in (xyz, velocity, weights, Bcells)) or np.any(weights <= 0):
        raise ValueError("Nonfinite data or nonpositive particle weights")
    if len(weights) and np.linalg.norm(velocity, axis=1).max() > .05*C:
        raise ValueError("Nonrelativistic ion moments require v < 0.05 c")
    ids = grid.ids(xyz)
    keep = np.isin(ids, region["cell_ids"])
    ids, v, w = ids[keep], velocity[keep], weights[keep]
    selected = np.array(region["cell_ids"])
    volume = region["volume_m3"]
    W = np.bincount(ids, weights=w, minlength=n)
    W2 = np.bincount(ids, weights=w*w, minlength=n)
    counts = np.bincount(ids, minlength=n)
    mean = np.zeros((n, 3))
    for j in range(3):
        np.divide(np.bincount(ids, weights=w*v[:, j], minlength=n), W,
                  out=mean[:, j], where=W > 0)
    bnorm = np.linalg.norm(Bcells, axis=1)
    if np.any(bnorm[selected] < 1e-12):
        raise ValueError("Near-zero B in ROI: parallel/perpendicular direction undefined")
    b = np.zeros_like(Bcells)
    np.divide(Bcells, bnorm[:, None], out=b, where=bnorm[:, None] > 0)
    c = v-mean[ids]
    cpar = np.einsum("ij,ij->i", c, b[ids])
    cperp2 = np.maximum(0, np.einsum("ij,ij->i", c, c)-cpar*cpar)
    Ppar = MP*float(w @ (cpar*cpar))/volume
    Pperp = .5*MP*float(w @ cperp2)/volume
    total = float(w.sum())
    density = total/volume
    neff = np.zeros(n)
    np.divide(W*W, W2, out=neff, where=W2 > 0)
    poor = neff[selected] < 10
    result = dict(macro_count=len(w), physical_count=total, density_m3=density,
                  effective_count=total**2/float(w @ w) if len(w) else 0.,
                  cell_count=len(selected), occupied_cells=int(np.count_nonzero(W[selected])),
                  minimum_cell_markers=int(counts[selected].min()),
                  low_neff_cell_fraction=float(poor.mean()),
                  low_neff_weight_fraction=float(W[selected][poor].sum()/total) if total else None,
                  P_parallel_Pa=Ppar, P_perp_Pa=Pperp, P_scalar_Pa=(Ppar+2*Pperp)/3,
                  T_parallel_eV=Ppar/(density*QE) if total else None,
                  T_perp_eV=Pperp/(density*QE) if total else None,
                  T_scalar_eV=(Ppar+2*Pperp)/(3*density*QE) if total else None)
    # Absolute lab-frame velocities for PDFs; moments above are drift-subtracted.
    vpar = np.einsum("ij,ij->i", v, b[ids])
    vperp = np.sqrt(np.maximum(0, np.einsum("ij,ij->i", v, v)-vpar*vpar))
    return result, (vpar, vperp, w)


def histogram(vpar, vperp, weights, parallel_edges, perp_edges):
    H, _, _ = np.histogram2d(vpar, vperp, (parallel_edges, perp_edges), weights=weights)
    total = float(weights.sum())
    density = H / (total*np.diff(parallel_edges)[:, None]*np.diff(perp_edges)[None, :]) if total else H
    return density, float(H.sum()/total) if total else None


def fluid_moments(Te, Pe, rho, volume, gamma, initialized, n_floor):
    """Equal-volume cell averages; Pe/rho is distinct from volume-mean Te.

    Native nodal products are averaged to cells independently, so we do not
    enforce <Pe> = <rho><Te>. Iteration zero precedes pressure initialization.
    The regional internal energy is an inventory, NOT a conservation residual.
    """
    Te, Pe, rho = (np.asarray(a, float) for a in (Te, Pe, rho))
    if not all(np.isfinite(a).all() for a in (Te, Pe, rho)):
        raise ValueError("Nonfinite electron-fluid diagnostics")
    if initialized and (np.any(Te <= 0) or np.any(Pe < 0) or np.any(rho <= 0)):
        raise ValueError("Interior fluid diagnostics require positive T/density, nonnegative pressure")
    mean_rho = float(rho.mean())
    return dict(Te_eV=float(Te.mean()), Te_min_eV=float(Te.min()),
                Te_max_eV=float(Te.max()),
                Te_effective_eV=float(Pe.mean()/mean_rho) if initialized and mean_rho > 0 else None,
                Pe_Pa=float(Pe.mean()), ne_grid_m3=mean_rho/QE,
                electron_U_region_J=float(Pe.mean()*volume/(gamma-1))
                    if initialized and gamma > 1 else None,
                below_density_floor_fraction=float(np.mean(rho/QE < n_floor)) if initialized else None,
                initialized_pressure=bool(initialized))


def kinetic_energy(u, weights):
    """Stable relativistic particle K; u = momentum/(mass*c)."""
    u, weights = np.asarray(u, float), np.asarray(weights, float)
    if not np.isfinite(u).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Nonfinite momenta or nonpositive weights")
    u2 = np.einsum("ij,ij->i", u, u)
    return MP*C*C*float(weights @ (u2/(np.sqrt(1+u2)+1)))


def analytic_mirror(xyz):
    """Continuum curl of this fixture's A, used ONLY for pre-solver t=0 output.

    WarpX dumps iteration 0 before it adds the split external B. Later output
    contains the total simulated B and must never have this field added again.
    """
    x, y, z = np.asarray(xyz).T
    bmin, ratio, k = .001487675, 3.10951921, 2*np.pi/5
    avg, a = bmin*(ratio+1)/2, (ratio-1)/(ratio+1)
    return np.column_stack((.5*x*avg*a*k*np.sin(k*z),
                            .5*y*avg*a*k*np.sin(k*z),
                            avg*(1+a*np.cos(k*z))))


def field_xyz(series, name, iteration, coord=None):
    array, info = series.get_field(field=name, coord=coord, iteration=iteration)
    order = [list(info.axes.values()).index(axis) for axis in "xyz"]
    values = np.transpose(array, order)
    if values.shape != Grid().shape or not np.isfinite(values).all():
        raise ValueError(f"Unexpected field shape/nonfinite values for {name}: {values.shape}")
    for i, axis in enumerate("xyz"):
        expected = Grid().lo[i]+(np.arange(Grid().shape[i])+.5)*Grid().dx[i]
        if not np.allclose(getattr(info, axis), expected, rtol=0, atol=1e-10):
            raise ValueError("Fields must be at the expected simulation cell centers")
    return values.reshape(-1)


def analyze(run):
    from openpmd_viewer import OpenPMDTimeSeries
    manifest = json.loads((run/"manifest.json").read_text())
    execution = json.loads((run/"execution.json").read_text())
    if execution["returncode"] != 0:
        raise ValueError("Refusing to label an unsuccessful run as complete")
    grid = Grid()
    rois = regions(grid)
    if manifest["grid"] != grid.as_dict():
        raise ValueError("Run mesh does not match diagnostic grid")
    series = OpenPMDTimeSeries(str(run/"diags/diag/openpmd"), backend="h5py")
    if int(series.iterations[-1]) != manifest["steps"] or not np.isclose(
            series.t[-1], manifest["horizon_s"], rtol=1e-10, atol=1e-18):
        raise ValueError("Snapshots do not reach the requested final step/time")
    if set(series.avail_species) != {"protons", "background_ions"}:
        raise ValueError("Expected ONLY two ion particle populations, no kinetic electrons")
    rows, fluid, pdfs, audits = [], [], {}, []
    parallel_edges = np.linspace(-1e5, 1e5, 81)
    perp_edges = np.linspace(0, 1e5, 61)
    for iteration, t in zip(series.iterations, series.t):
        B = np.column_stack([field_xyz(series, "B", int(iteration), a) for a in "xyz"])
        reference = "saved total simulated B at cell centers"
        if iteration == 0 and np.max(np.abs(B)) < 1e-14:
            B = analytic_mirror(grid.centers())
            reference = "analytic initial applied B: t=0 dump precedes hybrid initialization"
        E = np.column_stack([field_xyz(series, "E", int(iteration), a) for a in "xyz"])
        Te = field_xyz(series, "Te", int(iteration))*KB/QE  # Native WarpX output is Kelvin.
        Pe = field_xyz(series, "Pe", int(iteration))
        rho = field_xyz(series, "rho", int(iteration))
        covered = field_xyz(series, "eb_covered", int(iteration))
        snap = dict(step=int(iteration), time_s=float(t), B_reference=reference,
                    max_E_V_m=float(np.linalg.norm(E, axis=1).max()),
                    max_B_T=float(np.linalg.norm(B, axis=1).max()), populations={})
        if manifest.get("electron_energy_equation") and iteration > 0:
            ue_names = manifest["electron_velocity_output_fields"]
            if len(ue_names) != 3:
                raise ValueError("Need three documented electron-fluid velocity outputs")
            ue = np.column_stack([field_xyz(series, name, int(iteration)) for name in ue_names])
            snap["electron_advection_cell_center_CFL_proxy"] = float(
                np.max(np.sum(np.abs(ue)/grid.dx, axis=1))*manifest["dt_s"])
            snap["CFL_caveat"] = "Averaged saved velocity, not a strict bound on internal nodal/marker CFL"
        for region, spec in rois.items():
            idx = spec["cell_ids"]
            if np.any(covered[idx] != 0):
                raise ValueError("ROI intersects embedded boundary: cell-volume treatment must change")
            # Nodal Pe and rho are averaged to centers separately; the relation is
            # exactly preserved for isothermal gamma=1, not general polytropes.
            closure_error = None
            if iteration > 0 and manifest["closure"] == "isothermal":
                closure_error = float(np.max(np.abs(Pe[idx]-rho[idx]*manifest["Te0_eV"]))/
                                      (manifest["n0_m3"]*QE*manifest["Te0_eV"]))
                if closure_error > 1e-8 or not np.allclose(Te[idx], manifest["Te0_eV"], atol=1e-8, rtol=0):
                    raise ValueError("Isothermal fluid diagnostic closure mismatch")
            fluid.append(dict(step=int(iteration), time_s=float(t), region=region,
                              **fluid_moments(Te[idx], Pe[idx], rho[idx], spec["volume_m3"],
                                              manifest["gamma"], iteration > 0, manifest["n_floor_m3"]),
                              closure_error=closure_error))
        for species in series.avail_species:
            arrays = series.get_particle(["x", "y", "z", "ux", "uy", "uz", "w"],
                                         species=species, iteration=int(iteration))
            xyz, u = np.column_stack(arrays[:3]), np.column_stack(arrays[3:6])
            velocity = C*u/np.sqrt(1+np.sum(u*u, axis=1))[:, None]
            w = np.asarray(arrays[6])
            snap["populations"][species] = dict(markers=len(w), physical_count=float(w.sum()),
                                               kinetic_energy_J=kinetic_energy(u, w))
            for region, spec in rois.items():
                result, velocities = moments(xyz, velocity, w, B, grid, spec)
                rows.append(dict(step=int(iteration), time_s=float(t), region=region, species=species, **result))
                H, fraction = histogram(*velocities, parallel_edges, perp_edges)
                rows[-1]["histogram_weight_fraction"] = fraction
                pdfs[f"{int(iteration)}_{region}_{species}"] = H
        audits.append(snap)
    for name, data in (("ion_history", rows), ("electron_fluid_history", fluid)):
        with (run/f"{name}.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    np.savez_compressed(run/"velocity_distributions.npz", parallel_edges_m_s=parallel_edges,
                        perp_edges_m_s=perp_edges, **pdfs)
    write_json(run/"sampling_cells.json", dict(grid=grid.as_dict(), regions=rois))
    write_json(run/"diagnostic_audit.json", dict(snapshots=audits,
        final_ions=rows[-4:], final_fluid=fluid[-2:],
        finite=True, no_kinetic_electrons=True,
        electron_closure=manifest["closure"],
        caveat="Startup diagnostics, not beam heating or steady state. Energy inventories are not a closed budget for this open/source-driven mirror."))
    plot(run, rows, fluid, pdfs, parallel_edges, perp_edges, rois)
    print(json.dumps(dict(final_ions=rows[-4:], final_fluid=fluid[-2:]), indent=2))


def plot(run, rows, fluid, pdfs, par_edges, perp_edges, rois):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for i, region in enumerate(rois):
        selected = [r for r in rows if r["region"] == region and r["species"] == "background_ions"]
        ts = [r["time_s"]*1e9 for r in selected]
        for name, label in (("parallel", "parallel"), ("perp", "perpendicular"), ("scalar", "scalar")):
            axes[i, 0].plot(ts, [r[f"T_{name}_eV"] for r in selected], label=label)
            axes[i, 1].plot(ts, [r[f"P_{name}_Pa"]*1e3 for r in selected], label=label)
        f = [r for r in fluid if r["region"] == region and r["initialized_pressure"]]
        axes[i, 2].plot([r["time_s"]*1e9 for r in f], [r["Te_eV"] for r in f], color="C3")
        # Use physical scales, not auto-offset magnification of startup noise.
        axes[i, 0].set_ylim(0, max(1.2, 1.1*max(r[f"T_{key}_eV"] for r in selected
                                              for key in ("parallel", "perp", "scalar"))))
        axes[i, 1].set_ylim(0, max(1.9, 1.1e3*max(r[f"P_{key}_Pa"] for r in selected
                                               for key in ("parallel", "perp", "scalar"))))
        te_values = [r["Te_eV"] for r in f]
        axes[i, 2].set_ylim(min(9, min(te_values)*.9), max(11, max(te_values)*1.1))
        for j, ylabel in enumerate(("Background ion T (eV)", "Background ion P (mPa)", "Fluid Te (eV): see run closure")):
            axes[i, j].set(xlabel="Time (ns)", ylabel=ylabel, title=region.replace("_", " "))
            axes[i, j].grid(alpha=.2)
        axes[i, 0].legend()
    fig.suptitle("True hybrid-PIC | analytic periodic mirror TEST\nCell-local drift removed; startup changes are NOT beam-heating evidence")
    fig.savefig(run/"ion_fluid_history.png", dpi=170)
    plt.close(fig)
    last = max(r["step"] for r in rows)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    for i, region in enumerate(rois):
        for j, species in enumerate(("background_ions", "protons")):
            H = pdfs[f"{last}_{region}_{species}"]
            axes[i, j].set(title=f"{region} | {species}", xlabel="v parallel (km/s)", ylabel="v perpendicular (km/s)")
            if not H.any():
                axes[i, j].set(xlim=(par_edges[0]/1e3, par_edges[-1]/1e3),
                               ylim=(perp_edges[0]/1e3, perp_edges[-1]/1e3), facecolor="#f1f3f5")
                axes[i, j].text(.5, .5, "No particles in this region", ha="center", transform=axes[i, j].transAxes)
                continue  # No fabricated density scale for an undefined PDF.
            vmax = max(pdfs[f"{last}_{r}_{species}"].max()*1e6 for r in rois)
            im = axes[i, j].pcolormesh(par_edges/1e3, perp_edges/1e3, H.T*1e6,
                                      shading="auto", vmin=0, vmax=vmax)
            fig.colorbar(im, ax=axes[i, j], label="Probability density / (km/s)^2")
    fig.suptitle("Final ion velocity distributions | local total B; lab-frame velocities\nEach population separately physical-weight normalized; radial v-perp PDF")
    fig.savefig(run/"velocity_distributions.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for i, region in enumerate(rois):
        for species in ("background_ions", "protons"):
            selected = [r for r in rows if r["region"] == region and r["species"] == species]
            ts = [r["time_s"]*1e9 for r in selected]
            for j, key in enumerate(("density_m3", "macro_count", "physical_count")):
                axes[i, j].plot(ts, [r[key] for r in selected], label=species)
        for j, label in enumerate(("Density (1/m^3)", "Computational markers", "Represented physical particles")):
            axes[i, j].set(xlabel="Time (ns)", ylabel=label, title=region.replace("_", " "), yscale="symlog")
            axes[i, j].grid(alpha=.2)
        axes[i, 0].legend()
    fig.suptitle("Counts are not interchangeable | symlog axes include zero\nInherited beam is extremely weak; many markers can represent <<1 physical particle")
    fig.savefig(run/"density_and_counts.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for region, spec in rois.items():
        ijk = np.array(spec["ijk"])
        centers = np.array(Grid().lo)+(ijk+.5)*Grid().dx
        for ax, a, b in ((axes[0], 0, 1), (axes[1], 2, 0)):
            unique = np.unique(centers[:, [a, b]], axis=0)
            ax.scatter(unique[:, 0], unique[:, 1], marker="s", s=34, alpha=.65, label=region)
    axes[0].set(xlabel="x (m)", ylabel="y (m)", title="Selected cells: transverse projection", aspect="equal")
    axes[1].set(xlabel="z (m)", ylabel="x (m)", title="Selected cells: axial projection")
    for ax in axes:
        ax.legend(); ax.grid(alpha=.2)
    fig.suptitle("Exact cell-center selections; full cell volumes, zero-based indices saved in sampling_cells.json")
    fig.savefig(run/"sampling_cells.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    analyze(parser.parse_args().run)
