"""Input-generation guards; no WarpX execution or output-file creation."""
import math
import unittest

from setup_case import build_input


def parameters(deck):
    result = {}
    for line in deck.splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            result[key.strip()] = value.strip()
    return result


class SetupCaseTests(unittest.TestCase):
    def test_legacy_default_stays_isothermal(self):
        text, manifest = build_input()
        p = parameters(text)
        self.assertEqual(manifest["closure"], "isothermal")
        self.assertFalse(manifest["electron_energy_equation"])
        self.assertFalse(manifest["electron_temperature_evolves"])
        self.assertFalse(manifest["adaptive_fields"])
        self.assertEqual(p["hybrid_pic_model.solve_electron_energy_equation"], "0")
        self.assertEqual(p["diag.intervals"], "5")
        self.assertNotIn("hybrid_pic_model.use_rkf45", p)
        self.assertNotIn("diag.additional_fields_to_plot", p)

    def test_polytropic_is_not_transport(self):
        text, manifest = build_input(closure="polytropic")
        self.assertEqual(manifest["gamma"], 5/3)
        self.assertTrue(manifest["electron_temperature_evolves"])
        self.assertFalse(manifest["electron_energy_equation"])
        self.assertEqual(parameters(text)["hybrid_pic_model.solve_electron_energy_equation"], "0")

    def test_energy_transport_configuration(self):
        text, manifest = build_input(closure="energy_transport")
        p = parameters(text)
        self.assertEqual(manifest["gamma"], 5/3)
        self.assertTrue(manifest["electron_energy_equation"])
        self.assertTrue(manifest["electron_temperature_evolves"])
        self.assertFalse(manifest["electron_heating"])
        self.assertEqual(p["hybrid_pic_model.solve_electron_energy_equation"], "1")
        self.assertEqual(p["hybrid_pic_model.include_joule_heating"], "0")
        self.assertEqual(p["hybrid_pic_model.plasma_resistivity(rho,J,t)"], '"0.0"')
        self.assertFalse(any("relaxation_rate" in key for key in p))
        self.assertEqual(p["diag.additional_fields_to_plot"].split(), [
            '"hybrid_electron_velocity_fp[dir=x]"', '"hybrid_electron_velocity_fp[dir=y]"',
            '"hybrid_electron_velocity_fp[dir=z]"'])
        self.assertEqual(manifest["electron_velocity_output_fields"], [
            "hybrid_electron_velocity_fp[dir=x]", "hybrid_electron_velocity_fp[dir=y]",
            "hybrid_electron_velocity_fp[dir=z]"])
        self.assertIn("Te is K", text)
        self.assertTrue(any("no conduction" in item for item in manifest["limits"]))

    def test_energy_does_not_change_source_or_geometry(self):
        baseline, before = build_input()
        upgraded, after = build_input(closure="energy_transport")
        a, b = parameters(baseline), parameters(upgraded)
        prefixes = ("protons.", "background_ions.", "my_constants.", "geometry.",
                    "boundary.", "external_vector_potential.", "amr.")
        for key in a:
            if key.startswith(prefixes):
                self.assertEqual(a[key], b[key], key)
        self.assertEqual(a["warpx.eb_implicit_function"], b["warpx.eb_implicit_function"])
        self.assertEqual(before["beam_rate_physical_s"], after["beam_rate_physical_s"])
        self.assertEqual(before["source_sha256"], after["source_sha256"])
        self.assertEqual(before["regions"], after["regions"])

    def test_adaptive_controls_and_timestep_limit(self):
        text, manifest = build_input(dt=5e-9, adaptive_fields=True)
        p = parameters(text)
        self.assertEqual(p["hybrid_pic_model.use_rkf45"], "1")
        self.assertAlmostEqual(float(p["hybrid_pic_model.substep_rtol"]), 1e-5)
        self.assertAlmostEqual(float(p["hybrid_pic_model.substep_atol"]), 1e-12)
        self.assertTrue(manifest["adaptive_fields"])
        self.assertEqual(manifest["field_substep_rtol"], 1e-5)
        self.assertEqual(manifest["field_substep_atol"], 1e-12)
        self.assertEqual(manifest["dt_s"], 5e-9)
        for dt in (1.01e-9, 5e-9):
            with self.subTest(dt=dt), self.assertRaises(ValueError):
                build_input(dt=dt)
        with self.assertRaises(ValueError):
            build_input(dt=5.01e-9, adaptive_fields=True)

    def test_invalid_timesteps(self):
        for adaptive in (False, True):
            for dt in (0, -1e-9, math.nan, math.inf):
                with self.subTest(dt=dt, adaptive=adaptive), self.assertRaises(ValueError):
                    build_input(dt=dt, adaptive_fields=adaptive)

    def test_explicit_diagnostic_cadence(self):
        text, manifest = build_input(steps=41, diagnostic_every=7)
        self.assertEqual(parameters(text)["diag.intervals"], "7")
        self.assertEqual(manifest["diagnostic_every_steps"], 7)
        default, _ = build_input(steps=41)
        self.assertEqual(parameters(default)["diag.intervals"], "10")

    def test_invalid_cadence(self):
        for cadence in (0, -1, 1.5, True, "2"):
            with self.subTest(cadence=cadence), self.assertRaises(ValueError):
                build_input(diagnostic_every=cadence)

    def test_ion_collisions_do_not_enable_electron_exchange(self):
        text, manifest = build_input(closure="energy_transport", collisions=True)
        p = parameters(text)
        self.assertTrue(manifest["ion_collisions"])
        self.assertFalse(manifest["electron_heating"])
        self.assertEqual(p["collisions.collision_names"], "beam_ion ion_ion")
        self.assertFalse(any("relaxation_rate" in key for key in p))

    def test_unknown_closure_rejected(self):
        with self.assertRaises(ValueError):
            build_input(closure="kinetic_electrons")


if __name__ == "__main__":
    unittest.main()
