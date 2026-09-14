"""Input-composition and scale checks; these do not run WarpX or validate physics.

Run in the existing WSL environment with ``python -m unittest -v test_workflow``.
Generated inputs use temporary directories; all source assets remain read-only.
"""
from __future__ import annotations

import ast
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import workflow


def deck_values(text: str) -> dict[str, str]:
    """Reject duplicate active keys rather than silently taking the last value."""
    values = {}
    for line in text.splitlines():
        active = line.split("#", 1)[0].strip()
        if not active:
            continue
        key, value = active.split("=", 1)
        key, value = key.strip(), value.strip()
        if key in values:
            raise AssertionError(f"Duplicate active input key: {key}")
        values[key] = value
    return values


def evaluate_mask(expression: str, *, x: float, y: float, z: float) -> float:
    """Evaluate only elementary arithmetic/comparisons, never arbitrary code."""
    tree = ast.parse(expression.strip('"'), mode="eval")
    allowed = (ast.Expression, ast.BinOp, ast.Add, ast.Sub, ast.Mult, ast.Div,
               ast.Compare, ast.Lt, ast.Gt, ast.LtE, ast.GtE, ast.Name,
               ast.Load, ast.Constant, ast.UnaryOp, ast.USub)
    for node in ast.walk(tree):
        if not isinstance(node, allowed):
            raise AssertionError(f"Unsupported density expression node: {type(node).__name__}")
        if isinstance(node, ast.Name) and node.id not in {"x", "y", "z", "n_bg"}:
            raise AssertionError(f"Unknown density symbol: {node.id}")
    return float(eval(compile(tree, "<density-mask>", "eval"),
                      {"__builtins__": {}}, {"x": x, "y": y, "z": z, "n_bg": 1.0}))


@unittest.skipUnless(workflow.FIELD_DEFAULT.is_file(), "Local mirror HDF5 asset is not available")
class GeneratedModels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="mirror-background-test-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.before = {path: workflow.digest(path)
                      for path in (workflow.LARS, workflow.STL, workflow.FIELD_DEFAULT)}
        cls.original = deck_values(workflow.LARS.read_text(encoding="utf-8"))
        cls.generated = {}
        for case in workflow.MODELS:
            run = Path(cls.temporary.name) / case
            run.mkdir()
            metadata = workflow.generate(case, run, workflow.FIELD_DEFAULT, 20, 2, "analytic-bore")
            cls.generated[case] = (deck_values((run / "inputs.txt").read_text()), metadata)

    def test_source_assets_unchanged_and_live_source_retained(self):
        for path, before in self.before.items():
            with self.subTest(path=path):
                self.assertEqual(workflow.digest(path), before)
        retained = {key: value for key, value in self.original.items()
                    if key.startswith("protons.") and key not in
                    {"protons.do_not_deposit", "protons.initialize_self_fields"}}
        self.assertEqual(retained["protons.injection_style"], "NFluxPerCell")
        self.assertNotIn("protons.npart", retained)  # Commented Gaussian source stays inactive.
        for case, (deck, metadata) in self.generated.items():
            with self.subTest(case=case):
                for key, value in retained.items():
                    self.assertEqual(deck[key], value, key)
                for key in ("my_constants.y0", "my_constants.z0",
                            "my_constants.sigma_y", "my_constants.sigma_z"):
                    self.assertEqual(deck[key], self.original[key])
                self.assertEqual(metadata["source_sha256"], self.before[workflow.LARS])

    def test_explicit_standard_periodicity_and_output_species(self):
        for case, (deck, metadata) in self.generated.items():
            with self.subTest(case=case):
                for suffix in ("lo", "hi"):
                    self.assertEqual(deck[f"boundary.particle_{suffix}"], "absorbing absorbing periodic")
                    self.assertEqual(deck[f"boundary.field_{suffix}"], "pec pec periodic")
                self.assertFalse(any("Periodic_ReflectParVel" in value for value in deck.values()))
                self.assertEqual(deck["diag.species"].split(), metadata["species"])
                for species in metadata["species"]:
                    self.assertIn(f"diag.{species}.variables", deck)

    def test_kinetic_real_species_push_deposition_and_no_double_counting(self):
        for case in ("kinetic", "kinetic_half_dt"):
            deck, metadata = self.generated[case]
            with self.subTest(case=case):
                self.assertEqual(deck["particles.species_names"].split(),
                                 ["protons", "background_ions", "electrons"])
                self.assertEqual(deck["electrons.species_type"], "electron")
                self.assertEqual(deck["background_ions.species_type"], "proton")
                self.assertNotIn("electrons.mass", deck)
                self.assertEqual(metadata["electron_mass_kg"], workflow.ME)
                for species in metadata["species"]:
                    self.assertEqual(deck[f"{species}.do_not_deposit"], "0")
                    # WarpX defaults to pushing the inherited proton source.
                    self.assertEqual(deck.get(f"{species}.do_not_push", "0"), "0")
                self.assertTrue(metadata["self_consistent_E"])
                self.assertFalse(metadata["analytic_electron_drag"])
                names = deck["collisions.collision_names"].split()
                self.assertEqual(len(names), len(set(names)))
                pairs = []
                for name in names:
                    self.assertEqual(deck[f"{name}.type"], "pairwisecoulomb")
                    pairs.append(tuple(sorted(deck[f"{name}.species"].split())))
                self.assertEqual(len(pairs), len(set(pairs)))
                self.assertEqual(set(pairs), {
                    ("background_ions", "protons"), ("electrons", "protons"),
                    ("background_ions", "electrons"), ("background_ions", "background_ions"),
                    ("electrons", "electrons"), ("protons", "protons")})
                self.assertNotIn("background_stopping", deck.values())

    def test_reservoir_models_do_not_accidentally_enable_kinetic_electrons(self):
        for case in ("ion_background", "ion_coulomb", "hybrid"):
            deck, metadata = self.generated[case]
            with self.subTest(case=case):
                self.assertEqual(deck["particles.species_names"], "protons background_ions")
                self.assertEqual(deck["protons.do_not_deposit"], "1")
                self.assertEqual(deck["background_ions.do_not_deposit"], "1")
                self.assertEqual(deck["background_ions.do_not_push"], "1")
                self.assertFalse(metadata["self_consistent_E"])
                self.assertFalse(any(key.startswith("electrons.") for key in deck))
        control, _ = self.generated["ion_background"]
        ion, _ = self.generated["ion_coulomb"]
        hybrid, _ = self.generated["hybrid"]
        self.assertNotIn("collisions.collision_names", control)
        self.assertEqual(ion["collisions.collision_names"], "beam_ion")
        self.assertEqual(hybrid["collisions.collision_names"].split(), ["beam_ion", "electron_drag"])
        self.assertEqual(hybrid["electron_drag.type"], "background_stopping")
        self.assertEqual(hybrid["electron_drag.species"], "protons")
        self.assertEqual(hybrid["electron_drag.background_mass"], "m_e")
        self.assertEqual(hybrid["electron_drag.background_temperature"], "10*q_e/kb")

    def test_matching_background_masks_exclude_walls_and_periodic_seam(self):
        kinetic, _ = self.generated["kinetic"]
        hybrid, _ = self.generated["hybrid"]
        expressions = [kinetic["background_ions.density_function(x,y,z)"],
                       kinetic["electrons.density_function(x,y,z)"],
                       hybrid["background_ions.density_function(x,y,z)"],
                       hybrid["electron_drag.background_density(x,y,z,t)"]]
        self.assertEqual(len(set(expressions)), 1)
        for point, expected in (((0, 0, 2.5), 1), ((-0.2, 0, 1), 1),
                                ((0.749, 0, 1), 1), ((0.75, 0, 1), 0),
                                ((0.6, 0.6, 2.5), 0), ((0, 0, 0.25), 0),
                                ((0, 0, 4.75), 0), ((0, 0, 0), 0), ((0, 0, 5), 0)):
            with self.subTest(point=point):
                self.assertEqual(evaluate_mask(expressions[0], x=point[0], y=point[1], z=point[2]), expected)

    def test_electron_thermal_scale_and_combined_debye_resolution(self):
        deck, metadata = self.generated["kinetic"]
        electron_u = math.sqrt(10 * workflow.QE / workflow.ME) / workflow.C
        for axis in "xyz":
            self.assertAlmostEqual(float(deck[f"electrons.u{axis}_th"]) / electron_u, 1.0, places=13)
        # Inverse squared Debye lengths add for Ti=1 eV, Te=10 eV, equal densities.
        expected_lambda = math.sqrt(workflow.EPS0 / (
            metadata["density_m3"] * workflow.QE * (1 / 10 + 1 / 1)))
        self.assertAlmostEqual(metadata["lambda_D_combined_m"] / expected_lambda, 1.0, places=13)
        self.assertGreaterEqual(expected_lambda, metadata["cell_size_m"])
        self.assertLessEqual(metadata["electron_omega_p_dt"], 0.1)
        self.assertLessEqual(metadata["electron_omega_c_max_dt"], 0.1)

    def test_half_timestep_matches_physical_horizon_and_snapshot_cadence(self):
        full, a = self.generated["kinetic"]
        half, b = self.generated["kinetic_half_dt"]
        self.assertEqual(b["dt_s"], a["dt_s"] / 2)
        self.assertEqual(b["max_step"], a["max_step"] * 2)
        self.assertEqual(b["horizon_s"], a["horizon_s"])
        self.assertEqual(float(full["warpx.const_dt"]) * int(full["diag.intervals"]),
                         float(half["warpx.const_dt"]) * int(half["diag.intervals"]))
        for key in ("density_m3", "Ti_eV", "Te_eV", "ppc_each_dimension", "wall", "collision_names"):
            self.assertEqual(a[key], b[key], key)


class SourceChangeGuard(unittest.TestCase):
    def test_changed_source_style_requires_inspection(self):
        with tempfile.TemporaryDirectory(prefix="mirror-source-guard-") as directory:
            source = Path(directory) / "changed.txt"
            source.write_text("protons.injection_style = gaussian_beam\n", encoding="utf-8")
            with mock.patch.object(workflow, "LARS", source):
                with self.assertRaisesRegex(ValueError, "source changed"):
                    workflow.source_parameters()


if __name__ == "__main__":
    unittest.main()
