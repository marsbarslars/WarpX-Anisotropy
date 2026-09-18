import unittest

import numpy as np

from verify_energy import C, MP, QE, build_input, particle_budget, relative_change, thermal_reference


class EnergyVerificationTests(unittest.TestCase):
    def test_thermal_reference_initial_and_equilibrium(self):
        te, ti = thermal_reference([0., 1.], 10., 1., 1e5)
        np.testing.assert_allclose(te, [10., 5.5])
        np.testing.assert_allclose(ti, [1., 5.5])

    def test_exact_difference_rate_and_heat_capacity(self):
        t = np.array([0., 1e-7, 1e-6])
        te, ti = thermal_reference(t, 10., 1., 1e5)
        np.testing.assert_allclose(te-ti, 9*np.exp(-4e5*t))
        np.testing.assert_allclose(te+ti, 11.)

    def test_zero_rate_reference(self):
        te, ti = thermal_reference([0., 1.], 10., 1., 0.)
        np.testing.assert_allclose(te, 10.); np.testing.assert_allclose(ti, 1.)

    def test_thermal_deck_no_particles_electrons(self):
        deck, manifest = build_input()
        self.assertIn("particles.species_names = ions", deck)
        self.assertIn("solve_electron_energy_equation = 1", deck)
        self.assertIn("electron_ion_relaxation_rate(rho,Te,Ti,t)", deck)
        self.assertIn("diag.intervals = 1:1,10,100:100", deck)
        self.assertFalse(manifest["physical_collision_rate"])
        self.assertAlmostEqual(manifest["horizon_s"], 100e-9)

    def test_counterstream_has_equal_and_opposite_drift(self):
        deck, manifest = build_input("counterstream")
        parsed = dict(line.split(" = ", 1) for line in deck.splitlines() if " = " in line)
        self.assertAlmostEqual(float(parsed["plus.ux_m"]), -float(parsed["minus.ux_m"]))
        self.assertEqual(float(parsed["plus.density"]), 5e15)
        self.assertEqual(float(parsed["minus.density"]), 5e15)
        self.assertEqual(manifest["Te0_eV"], manifest["Ti0_eV"])

    def test_bad_parameters(self):
        for parameters in (dict(case="beam"), dict(steps=0), dict(dt=0), dict(dt=float("nan")),
                           dict(dt=1e-6), dict(nu=-1), dict(nu=1e9), dict(cells=7), dict(ppc=1)):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                build_input(**parameters)

    def test_particle_temperature_drift_and_energy(self):
        v = np.vstack((np.eye(3), -np.eye(3)))*3000.+[1e4, 0, 0]
        u = v/C/np.sqrt(1-np.sum((v/C)**2, axis=1))[:, None]
        result = particle_budget(u, np.ones(6)*2.)
        self.assertAlmostEqual(result["Ti_eV"], MP*3e6/QE)
        self.assertAlmostEqual(result["mean_vx_m_s"], 1e4)
        self.assertEqual(result["physical_count"], 12.)
        np.testing.assert_allclose(result["K_J"], result["bulk_K_J"]+1.5*12*QE*result["Ti_eV"], rtol=1e-7)

    def test_particle_weight_scaling(self):
        u = np.array([[1e-4, 0, 0], [-1e-4, 0, 0]])
        a, b = particle_budget(u, [1, 1]), particle_budget(u, [7, 7])
        self.assertAlmostEqual(a["Ti_eV"], b["Ti_eV"])
        np.testing.assert_allclose(b["K_J"], 7*a["K_J"], rtol=1e-14)

    def test_reject_bad_particle_data(self):
        for u, w in ((np.zeros((0, 3)), []), (np.zeros((1, 3)), [-1]),
                     (np.full((1, 3), np.nan), [1]), (np.ones((1, 3)), [1])):
            with self.assertRaises(ValueError): particle_budget(u, w)

    def test_signed_energy_residual(self):
        np.testing.assert_allclose(relative_change([10., 9., 11.]), [0., -.1, .1])
        with self.assertRaises(ValueError): relative_change([0., 1.])


if __name__ == "__main__":
    unittest.main()
