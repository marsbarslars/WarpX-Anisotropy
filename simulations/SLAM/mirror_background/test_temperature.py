"""Unit tests for temperature units, physical weights, drift, and probe inputs."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

import temperature
import workflow
from test_workflow import deck_values


class TemperatureMoments(unittest.TestCase):
    def test_known_isotropic_velocity_spread(self):
        speed = np.sqrt(3 * workflow.QE / workflow.MP)
        v = np.vstack((np.eye(3), -np.eye(3))) * speed
        result = temperature.thermal_moments(v, np.ones(6), workflow.MP)
        self.assertAlmostEqual(result["T_eV"], 1.0, places=13)
        self.assertAlmostEqual(result["mean_kinetic_energy_eV"], 1.5, places=13)
        self.assertEqual(result["bulk_energy_eV"], 0.0)

    def test_uniform_drift_is_not_heating(self):
        v = np.array([[1000., 2000., -4000.], [-1000., -2000., 4000.]])
        w = np.array([2., 5.])
        baseline = temperature.thermal_moments(v, w, workflow.MP)
        drifted = temperature.thermal_moments(v + [100000., 50000., -20000.], w, workflow.MP)
        self.assertAlmostEqual(drifted["T_eV"], baseline["T_eV"], places=12)
        self.assertGreater(drifted["bulk_energy_eV"], baseline["bulk_energy_eV"])

    def test_weights_and_macroparticle_splitting(self):
        v = np.array([[-3000., 0, 0], [9000., 0, 0]])
        weighted = temperature.thermal_moments(v, np.array([3., 1.]), workflow.MP)
        split = temperature.thermal_moments(np.repeat(v, [3, 1], axis=0), np.ones(4), workflow.MP)
        self.assertEqual(weighted["vx_mean_m_s"], 0)
        self.assertAlmostEqual(weighted["T_eV"], split["T_eV"], places=13)
        self.assertEqual(weighted["physical_weight"], split["physical_weight"])
        self.assertNotEqual(weighted["macro_count"], split["macro_count"])

    def test_empty_and_cold_region_are_distinct(self):
        empty = temperature.thermal_moments(np.empty((0, 3)), np.empty(0), workflow.MP)
        cold = temperature.thermal_moments(np.array([[1000., 0, 0]]), np.ones(1), workflow.MP)
        self.assertIsNone(empty["T_eV"])
        self.assertEqual(cold["T_eV"], 0)
        self.assertGreater(cold["bulk_energy_eV"], 0)

    def test_normalized_momentum_is_not_velocity(self):
        u = np.array([[0.01, 0.02, 0.03]])
        v = temperature.velocity_from_u(u)
        expected = workflow.C * u / np.sqrt(1.0014)
        np.testing.assert_allclose(v, expected, rtol=1e-15)

    def test_invalid_or_relativistic_particles_rejected(self):
        for w in (np.array([-1.]), np.array([0.]), np.array([np.nan])):
            with self.assertRaises(ValueError):
                temperature.thermal_moments(np.zeros((1, 3)), w, workflow.MP)
        with self.assertRaises(ValueError):
            temperature.thermal_moments(np.array([[workflow.C*0.1, 0, 0]]), np.ones(1), workflow.MP)

    def test_central_region_has_strict_bounds(self):
        probe = {"radius_m": 0.25, "z_min_m": 2.25, "z_max_m": 2.75}
        xyz = np.array([[0, 0, 2.5], [0.25, 0, 2.5], [0, 0, 2.25], [0, 0, 1]])
        np.testing.assert_array_equal(temperature.region_mask(*xyz.T, probe), [True, False, False, False])


@unittest.skipUnless(workflow.FIELD_DEFAULT.is_file(), "Local mirror field asset not available")
class TemperatureDeck(unittest.TestCase):
    def test_shareable_deck_matches_generated_kinetic_configuration(self):
        with tempfile.TemporaryDirectory(prefix="temperature-shareable-") as directory:
            run = Path(directory)
            workflow.generate("kinetic", run, workflow.FIELD_DEFAULT, 600, 2, "analytic-bore", 1)
            generated = deck_values((run / "inputs.txt").read_text())
            shared = deck_values((workflow.HERE / "inputs_kinetic_temperature_30ns.txt").read_text())
            generated["particles.read_fields_from_path"] = '"example-femm-3d.h5"'
            self.assertEqual(generated, shared)

    def test_probe_changes_diagnostics_not_physics(self):
        with tempfile.TemporaryDirectory(prefix="temperature-deck-") as directory:
            plain, probe = Path(directory)/"plain", Path(directory)/"probe"
            plain.mkdir()
            probe.mkdir()
            for case in workflow.MODELS:
                workflow.generate(case, plain, workflow.FIELD_DEFAULT, 20, 2, "analytic-bore")
                metadata = workflow.generate(case, probe, workflow.FIELD_DEFAULT, 20, 2, "analytic-bore", 1)
                a, b = (deck_values((folder/"inputs.txt").read_text()) for folder in (plain, probe))
                for key in a:
                    if key != "diagnostics.diags_names":
                        self.assertEqual(a[key], b[key], (case, key))
                self.assertEqual(b["diagnostics.diags_names"], "diag center")
                self.assertEqual(b["center.intervals"], "1")
                self.assertEqual(b["center.fields_to_plot"], "none")
                for species in metadata["species"]:
                    self.assertEqual(b[f"center.{species}.plot_filter_function(t,x,y,z,ux,uy,uz)"].strip('"'),
                                     metadata["temperature_probe"]["filter"])
                    self.assertNotIn(f"center.{species}.random_fraction", b)


if __name__ == "__main__":
    unittest.main()
