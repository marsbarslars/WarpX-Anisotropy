import unittest

import numpy as np

from diagnostics import C, Grid, MP, QE, analytic_mirror, fluid_moments, histogram, kinetic_energy, moments, regions
from setup_case import build_input


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.grid = Grid(lo=(0., 0., 0.), hi=(2., 1., 1.), shape=(2, 1, 1))
        self.B = np.array([[0., 0., 1.], [0., 0., 1.]])
        self.roi = dict(cell_ids=[0, 1], volume_m3=2.)
        self.xyz = np.array([[.5, .5, .5]]*6+[[1.5, .5, .5]]*6)
        self.v = np.tile(np.vstack((np.eye(3), -np.eye(3))), (2, 1))*3000.
        self.w = np.ones(12)*2.

    def calculate(self, v=None, w=None, B=None):
        return moments(self.xyz, self.v if v is None else v, self.w if w is None else w,
                       self.B if B is None else B, self.grid, self.roi)[0]

    def test_isotropy_and_pressure(self):
        r = self.calculate()
        self.assertAlmostEqual(r["T_parallel_eV"], MP*3e6/QE)
        self.assertAlmostEqual(r["T_perp_eV"], r["T_parallel_eV"])
        self.assertAlmostEqual(r["P_scalar_Pa"], r["density_m3"]*QE*r["T_scalar_eV"])
        self.assertEqual(r["density_m3"], 12.)

    def test_cell_local_not_region_flow(self):
        v = self.v.copy()
        v[:6] += (3e4, 2e4, -4e4)
        v[6:] -= (5e4, 3e4, 2e4)
        self.assertAlmostEqual(self.calculate(v=v)["T_scalar_eV"], self.calculate()["T_scalar_eV"])

    def test_local_field_rotation(self):
        v = np.zeros_like(self.v)
        v[:6, 0] = (-2., -1., 0., 0., 1., 2.)
        v[6:, 1] = (-2., -1., 0., 0., 1., 2.)
        r = self.calculate(v=v, B=np.array([[1., 0., 0.], [0., 1., 0.]]))
        self.assertGreater(r["T_parallel_eV"], 0)
        self.assertAlmostEqual(r["T_perp_eV"], 0)

    def test_weight_scaling(self):
        a, b = self.calculate(), self.calculate(w=17*self.w)
        self.assertAlmostEqual(a["T_scalar_eV"], b["T_scalar_eV"])
        self.assertAlmostEqual(17*a["P_parallel_Pa"], b["P_parallel_Pa"])
        self.assertAlmostEqual(a["effective_count"], b["effective_count"])

    def test_unequal_weights(self):
        w = self.w.copy(); w[0] *= 10
        r = self.calculate(w=w)
        self.assertAlmostEqual(r["effective_count"], w.sum()**2/(w@w))
        self.assertLess(r["effective_count"], len(w))

    def test_empty_not_cold(self):
        r, _ = moments(np.empty((0, 3)), np.empty((0, 3)), np.array([]), self.B, self.grid, self.roi)
        self.assertIsNone(r["T_scalar_eV"])
        self.assertEqual(r["P_scalar_Pa"], 0.)
        self.assertEqual(r["density_m3"], 0.)

    def test_single_marker_flag(self):
        r, _ = moments(self.xyz[:1], self.v[:1], self.w[:1], self.B, self.grid, self.roi)
        self.assertEqual(r["low_neff_weight_fraction"], 1.)

    def test_invalid_weight(self):
        w = self.w.copy(); w[0] = -1
        with self.assertRaises(ValueError): self.calculate(w=w)

    def test_nonfinite(self):
        v = self.v.copy(); v[0, 0] = np.nan
        with self.assertRaises(ValueError): self.calculate(v=v)

    def test_zero_B(self):
        with self.assertRaises(ValueError): self.calculate(B=np.zeros_like(self.B))

    def test_relativistic_rejected(self):
        with self.assertRaises(ValueError): self.calculate(v=np.ones_like(self.v)*C)

    def test_grid_boundaries(self):
        ids = self.grid.ids(np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.], [-.001, 0., 0.]]))
        np.testing.assert_array_equal(ids, (0, 1, -1, -1))

    def test_roi_cells(self):
        r = regions(Grid())
        self.assertEqual(r["center"]["count"], 416)
        self.assertEqual(r["near_injection"]["count"], 64)
        self.assertAlmostEqual(r["center"]["volume_m3"], .1015625)
        self.assertAlmostEqual(r["near_injection"]["volume_m3"], .015625)
        self.assertFalse(set(r["center"]["cell_ids"]) & set(r["near_injection"]["cell_ids"]))
        self.assertEqual(r["near_injection"]["bounds_m"], [[-.25, -.125, .875], [0., .125, 1.125]])

    def test_changed_mesh_rejected(self):
        with self.assertRaises(ValueError): regions(self.grid)

    def test_histogram_normalization_and_tail(self):
        par, perp = np.array([-2., 0., 2.]), np.array([1., 1., 1.])
        edges = np.array([-1., 1.]); edges2 = np.array([0., 2.])
        H, fraction = histogram(par, perp, np.ones(3), edges, edges2)
        self.assertAlmostEqual(fraction, 1/3)
        self.assertAlmostEqual(float(H.sum())*4, 1/3)

    def test_deck_is_true_hybrid(self):
        deck, manifest = build_input()
        self.assertIn("algo.maxwell_solver = hybrid", deck)
        self.assertIn("particles.species_names = protons background_ions", deck)
        self.assertNotIn("particles.B_ext_particle_init_style", deck)
        self.assertIn("hybrid_pic_model.add_external_fields = 1", deck)
        self.assertEqual(manifest["gamma"], 1.)
        self.assertFalse(manifest["electron_heating"])

    def test_polytropic_explicit(self):
        _, m = build_input(closure="polytropic")
        self.assertEqual(m["gamma"], 5/3)
        self.assertFalse(m["electron_heating"])

    def test_no_silent_timestep_increase(self):
        with self.assertRaises(ValueError): build_input(dt=1e-6)

    def test_applied_mirror_periodic_and_ratio(self):
        B = analytic_mirror(np.array([[0., 0., 0.], [0., 0., 2.5], [.2, .3, 0.], [.2, .3, 5.]]))
        self.assertAlmostEqual(B[0, 2]/B[1, 2], 3.10951921)
        np.testing.assert_allclose(B[2], B[3], atol=1e-15)

    def test_applied_mirror_divergence(self):
        p = np.array([[.2, .1, .9]])
        h = 1e-5
        divergence = 0.
        for axis in range(3):
            delta = np.zeros_like(p); delta[:, axis] = h
            divergence += (analytic_mirror(p+delta)[0, axis]-analytic_mirror(p-delta)[0, axis])/(2*h)
        self.assertLess(abs(divergence), 1e-11)

    def test_fluid_effective_temperature_not_volume_average(self):
        r = fluid_moments([1, 3], np.array([1, 9])*QE,
                          np.array([1, 3])*QE, 2., 5/3, True, .1)
        self.assertEqual(r["Te_eV"], 2.)
        self.assertAlmostEqual(r["Te_effective_eV"], 2.5)
        self.assertAlmostEqual(r["electron_U_region_J"], 15*QE)

    def test_uninitialized_pressure_not_cold_electrons(self):
        r = fluid_moments([10], [0], [0], 1, 5/3, False, 1e15)
        self.assertIsNone(r["Te_effective_eV"])
        self.assertIsNone(r["electron_U_region_J"])

    def test_negative_fluid_rejected(self):
        with self.assertRaises(ValueError):
            fluid_moments([10], [-1], [1], 1, 5/3, True, 1e15)

    def test_stable_kinetic_energy(self):
        u = np.array([[1e-8, 0., 0.]])
        self.assertAlmostEqual(kinetic_energy(u, [2])/(MP*C*C*1e-16), 1.)
        self.assertEqual(kinetic_energy(np.empty((0, 3)), []), 0.)


if __name__ == "__main__":
    unittest.main()
