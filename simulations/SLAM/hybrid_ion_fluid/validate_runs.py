"""Reanalyze the four named validation runs and compare identical-horizon dt tests.

This does not declare physical convergence of heating, losses, or the vessel.
All output is derived from completed runs, never synthesized time samples.
"""
from pathlib import Path
import json

import numpy as np

from diagnostics import analyze, write_json

HERE = Path(__file__).resolve().parent
NAMES = ("isothermal_20ns_20260918", "isothermal_halfdt_20ns_20260918",
         "isothermal_ioncoll_50ns_20260918", "polytropic_4ns_20260918")


def main():
    audits, summary = {}, {}
    for name in NAMES:
        path = HERE/"runs"/name
        analyze(path)
        audit = json.loads((path/"diagnostic_audit.json").read_text())
        manifest = json.loads((path/"manifest.json").read_text())
        execution = json.loads((path/"execution.json").read_text())
        audits[name] = audit
        background = [s["populations"]["background_ions"]["markers"] for s in audit["snapshots"]]
        assert min(background) == max(background) == 974400
        assert all(r["macro_count"] == 0 for r in audit["final_ions"]
                   if r["species"] == "protons" and r["region"] == "center")
        assert all(r["histogram_weight_fraction"] is None or
                   abs(r["histogram_weight_fraction"]-1) < 1e-10 for r in audit["final_ions"])
        summary[name] = dict(execution=execution, horizon_s=manifest["horizon_s"],
                             dt_s=manifest["dt_s"], ion_collisions=manifest["ion_collisions"],
                             closure=manifest["closure"], initial_background_markers=background[0],
                             final_ions=audit["final_ions"], final_fluid=audit["final_fluid"])
    a, b = (audits[n] for n in NAMES[:2])
    comparisons = []
    for region in ("center", "near_injection"):
        aa, bb = [next(r for r in result["final_ions"] if r["region"] == region and
                     r["species"] == "background_ions") for result in (a, b)]
        differences = {key: float(abs(aa[key]-bb[key])/abs(aa[key])) for key in
                       ("T_parallel_eV", "T_perp_eV", "T_scalar_eV", "P_scalar_Pa", "density_m3")}
        assert max(differences.values()) < 1e-4
        comparisons.append(dict(region=region, relative_differences=differences))
    write_json(HERE/"runs"/"validation_summary.json", dict(runs=summary,
               half_dt_comparison=comparisons, operational_checks_passed=True,
               caveat="Only short-startup background moments checked; not production convergence"))
    print("All completed-run, count, histogram and short-horizon background timestep checks passed.")
    print(json.dumps(comparisons, indent=2))


if __name__ == "__main__":
    main()
