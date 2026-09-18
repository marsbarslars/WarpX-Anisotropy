import copy
import csv
from pathlib import Path
import tempfile
import unittest

import numpy as np

from presentation_plots import QE, compare_background, load_pdf_payload, read_rows, selection, time_axis


def fixture(times=(0., 1e-8, 2e-8), dt=1e-9):
    manifest = dict(model="hybrid", geometry="periodic mirror", closure="energy_transport",
                    source_sha256="abc", grid={"shape": [2, 2, 2]}, n0_m3=1e16,
                    n_floor_m3=1e15, Te0_eV=10., Ti0_eV=1., ppc_per_axis=3,
                    field_substeps=64, ion_collisions=True, dt_s=dt,
                    steps=round(2e-8/dt), horizon_s=2e-8)
    ions = []
    for region in ("center", "near_injection"):
        for time in times:
            ions.append(dict(time_s=time, step=round(time/dt), region=region,
                             species="background_ions", T_parallel_eV=1., T_perp_eV=1.,
                             T_scalar_eV=1., P_parallel_Pa=.001, P_perp_Pa=.001,
                             P_scalar_Pa=.001, density_m3=1e16))
    fluid = [dict(time_s=time, step=round(time/dt), region=region, initialized_pressure=time > 0,
                  Te_eV=10., Te_effective_eV=10. if time > 0 else None,
                  Pe_Pa=1e16*QE*10 if time > 0 else 0., ne_grid_m3=1e16)
             for region in ("center", "near_injection") for time in times]
    return {"manifest": manifest, "ions": ions, "fluid": fluid}


class PresentationTests(unittest.TestCase):
    def test_zero_misfit(self):
        a, b = fixture(), fixture(dt=.5e-9)
        result = compare_background(a, b)
        self.assertEqual(result["max_abs_scaled_difference"], 0)
        self.assertFalse(result["interpolation"])
        self.assertFalse(result["reference_is_truth"])
        self.assertEqual(result["max_abs_fluid_scaled_difference"], 0)
        self.assertTrue(all(r["time_s"] > 0 for r in result["fluid_records"]))

    def test_fluid_normalization_and_separate_maximum(self):
        a, b = fixture(), fixture(dt=.5e-9)
        a["fluid"][-1]["Te_eV"] = 11.
        a["fluid"][-1]["Te_effective_eV"] = 12.
        a["fluid"][-1]["Pe_Pa"] *= 1.3
        a["fluid"][-1]["ne_grid_m3"] *= 1.4
        result = compare_background(a, b)
        last = {r["quantity"]: r["signed_scaled_difference"] for r in result["fluid_records"]
                if r["region"] == "near_injection" and r["time_s"] == 2e-8}
        self.assertAlmostEqual(last["Te_eV"], .1)
        self.assertAlmostEqual(last["Te_effective_eV"], .2)
        self.assertAlmostEqual(last["Pe_Pa"], .3)
        self.assertAlmostEqual(last["ne_grid_m3"], .4)
        self.assertEqual(result["max_abs_scaled_difference"], 0)
        self.assertAlmostEqual(result["max_abs_fluid_scaled_difference"], .4)

    def test_optional_effective_fluid_temperature(self):
        a, b = fixture(), fixture()
        for row in a["fluid"]:
            del row["Te_effective_eV"]
        result = compare_background(a, b)
        self.assertFalse(result["fluid_sampling"]["center"]["effective_temperature_compared"])
        self.assertNotIn("Te_effective_eV", {r["quantity"] for r in result["fluid_records"]})

    def test_uninitialized_fluid_endpoint_rejected(self):
        a, b = fixture(), fixture()
        b["fluid"][-1]["initialized_pressure"] = False
        with self.assertRaises(ValueError): compare_background(a, b)

    def test_signed_initial_scale_not_tiny_change(self):
        a, b = fixture(), fixture(dt=.5e-9)
        a["ions"][-1]["T_parallel_eV"] = 1.02
        a["ions"][-1]["P_scalar_Pa"] = .00099
        a["ions"][-1]["density_m3"] = 1.01e16
        records = compare_background(a, b)["records"]
        final = {r["quantity"]: r for r in records if r["region"] == "near_injection" and r["time_s"] == 2e-8}
        self.assertAlmostEqual(final["T_parallel_eV"]["signed_scaled_difference"], .02)
        self.assertAlmostEqual(final["P_scalar_Pa"]["signed_scaled_difference"], -.01)
        self.assertAlmostEqual(final["density_m3"]["signed_scaled_difference"], .01)

    def test_no_interpolation(self):
        a, b = fixture(), fixture(times=(0., .7e-8, 2e-8), dt=.5e-9)
        result = compare_background(a, b)
        self.assertEqual(result["sampling"]["center"]["matched_times_s"], [0., 2e-8])
        self.assertEqual(result["sampling"]["center"]["unmatched_run_samples"], 1)
        self.assertEqual(result["fluid_sampling"]["center"]["matched_times_s"], [2e-8])

    def test_roundoff_time_match(self):
        a, b = fixture(), fixture(times=(0., 1e-8+1e-23, 2e-8+2e-23))
        self.assertEqual(len(compare_background(a, b)["sampling"]["center"]["matched_times_s"]), 3)

    def test_diagnostic_cadence_is_not_a_physics_change(self):
        for key in ("diagnostic_every", "diagnostic_every_steps", "diagnostics_every", "diagnostics_every_steps"):
            with self.subTest(key=key):
                a, b = fixture(), fixture(dt=.5e-9)
                a["manifest"][key] = 10
                b["manifest"][key] = 20
                self.assertEqual(compare_background(a, b)["max_abs_scaled_difference"], 0)

    def test_scenario_changes_rejected(self):
        for key, value in (("closure", "isothermal"), ("source_sha256", "different"),
                           ("geometry", "SLAM"), ("n0_m3", 1e17),
                           ("ion_collisions", False), ("ppc_per_axis", 4),
                           ("electron_relaxation_rate", 2.)):
            with self.subTest(key=key):
                a, b = fixture(), fixture()
                b["manifest"][key] = value
                with self.assertRaises(ValueError): compare_background(a, b)

    def test_missing_metadata_rejected(self):
        a, b = fixture(), fixture()
        del b["manifest"]["source_sha256"]
        with self.assertRaises(ValueError): compare_background(a, b)

    def test_mismatched_horizon_rejected(self):
        a, b = fixture(), fixture()
        b["manifest"]["horizon_s"] = 3e-8
        with self.assertRaises(ValueError): compare_background(a, b)

    def test_missing_endpoint_rejected(self):
        with self.assertRaises(ValueError): compare_background(fixture(times=(0., 1e-8)), fixture())

    def test_undefined_or_zero_scale_rejected(self):
        a, b = fixture(), fixture()
        a["ions"][0]["T_scalar_eV"] = None
        with self.assertRaises(ValueError): compare_background(a, b)
        a, b = fixture(), fixture()
        b["ions"][0]["P_scalar_Pa"] = 0
        with self.assertRaises(ValueError): compare_background(a, b)

    def test_uninitialized_pressure_omitted(self):
        rows = [{"region": "center", "step": 0, "time_s": 0, "initialized_pressure": False},
                {"region": "center", "step": 1, "time_s": 1e-9, "initialized_pressure": True}]
        self.assertEqual(selection(rows, "center", initialized=True), rows[1:])

    def test_duplicate_times_rejected(self):
        rows = fixture()["ions"]
        rows.append(copy.deepcopy(rows[0]))
        with self.assertRaises(ValueError): selection(rows, "center", "background_ions")

    def test_csv_blank_is_undefined_not_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"data.csv"
            with path.open("w", newline="") as stream:
                writer = csv.writer(stream)
                writer.writerow(("region", "species", "Te_eV", "initialized_pressure"))
                writer.writerow(("center", "protons", "", "False"))
            row = read_rows(path)[0]
            self.assertIsNone(row["Te_eV"])
            self.assertIs(row["initialized_pressure"], False)

    def test_units_follow_horizon(self):
        self.assertEqual(time_axis(5e-8), (1e9, "ns"))
        self.assertEqual(time_axis(1e-5), (1e6, "microseconds"))

    def test_missing_legacy_pdf_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertIsNone(load_pdf_payload({"path": Path(directory)}))

    def test_pdf_probability_preserved_not_renormalized(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            # One bin with area 4 (m/s)^2 contains half the population weight.
            np.savez(path/"velocity_distributions.npz", parallel_edges_m_s=[-1., 1.],
                     perp_edges_m_s=[0., 2.],
                     **{"1_center_background_ions": np.array([[.125]]),
                        "1_near_injection_protons": np.array([[0.]])})
            ions = [dict(region="center", species="background_ions", time_s=1e-9, step=1,
                         macro_count=10, histogram_weight_fraction=.5),
                    dict(region="near_injection", species="protons", time_s=1e-9, step=1,
                         macro_count=0, histogram_weight_fraction=None),
                    dict(region="center", species="protons", time_s=1e-9, step=1,
                         macro_count=0, histogram_weight_fraction=None)]
            payload = load_pdf_payload({"path": path, "ions": ions})
            self.assertEqual(payload["panels"][0]["captured_weight_fraction"], .5)
            self.assertEqual(payload["panels"][0]["H"][0, 0], .125)
            self.assertEqual(payload["central_beam_macro_count"], 0)
            ions[0]["histogram_weight_fraction"] = .9
            with self.assertRaises(ValueError): load_pdf_payload({"path": path, "ions": ions})


if __name__ == "__main__":
    unittest.main()
