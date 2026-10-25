"""Small numerical and endpoint regressions; no competition-scale scan."""
import unittest
import numpy as np

from space_dust import historical
from space_dust.integrate import integrate
from space_dust.model import initial_state, invariant_residual, mass_from_radius, positive_rhs, radius_from_mass


class ModelTests(unittest.TestCase):
    def test_mass_units_and_radius_round_trip(self):
        mass = float(mass_from_radius(1e-4))
        self.assertAlmostEqual(mass * 1e6, 0.012566370614359173)
        self.assertAlmostEqual(float(radius_from_mass(mass)), 1e-4)

    def test_mass_and_radius_percentages_differ(self):
        result = integrate([100, 1, 1], t_max=3, dt=0.1,
                           rhs=lambda t, y: np.array([0, 0, -0.25]), mass_fraction_cutoff=0.5)
        summary = result.summary()
        self.assertEqual(result.status, "mass_cutoff")
        self.assertAlmostEqual(summary["end_mass_percent"], 50)
        self.assertAlmostEqual(summary["end_radius_percent"], 100 * 0.5 ** (1 / 3))
        self.assertIsNone(summary["ground_mass_percent"])

    def test_guarded_and_historical_rk4_agree_in_valid_short_case(self):
        initial = initial_state(0.1, 12)
        old = historical.solve_rk4_manual(historical.particle_model, initial, [0, 0.01], dt=0.001)
        current = integrate(initial, t_max=0.01, dt=0.001)
        np.testing.assert_allclose(current.states[-1], old, rtol=2e-13, atol=1e-12)
        self.assertEqual(current.rejected_steps, 0)

    def test_time_limit_has_no_ground_mass(self):
        result = integrate(initial_state(0.1, 12), t_max=0.01, dt=0.003)
        self.assertEqual(result.status, "time_limit")
        self.assertEqual(result.time_s[-1], 0.01)
        self.assertEqual(len(result.time_s), len(result.states))
        self.assertGreater(result.states[-1, 0], 0)
        self.assertIsNone(result.summary()["ground_mass_percent"])

    def test_ground_event_is_interpolated_and_arrays_match(self):
        result = integrate([1, 2, 1], t_max=2, dt=1, rhs=lambda t, y: np.array([-2, 0, 0]))
        self.assertEqual(result.status, "ground_reached")
        self.assertAlmostEqual(result.time_s[-1], 0.5)
        self.assertEqual(result.states[-1, 0], 0)
        self.assertEqual(len(result.time_s), len(result.states))
        self.assertEqual(result.summary()["ground_mass_percent"], 100)

    def test_nonfinite_rhs_returns_numerical_failure(self):
        result = integrate([100, 1, 1], t_max=1, dt=0.1,
                           rhs=lambda t, y: np.array([np.nan, 0, 0]), max_halvings=1)
        self.assertEqual(result.status, "numerical_failure")
        self.assertTrue(np.isfinite(result.states).all())

    def test_negative_velocity_is_never_accepted(self):
        result = integrate([100, 1, 1], t_max=3, dt=2,
                           rhs=lambda t, y: np.array([-1, -1, 0]), max_halvings=3)
        self.assertEqual(result.status, "numerical_failure")
        self.assertGreater(result.rejected_steps, 0)
        self.assertTrue((result.states[:, 1] >= 0).all())

    def test_small_particle_guard_is_bounded(self):
        result = integrate(initial_state(1e-6, 6), t_max=0.01, dt=0.01, max_steps=3)
        self.assertIn(result.status, {"step_limit", "numerical_failure", "time_limit"})
        self.assertGreater(result.rejected_steps, 0)
        self.assertTrue(np.isfinite(result.states).all())
        self.assertTrue((result.states[:, 1:] >= 0).all())

    def test_exact_identity_is_tangent_to_original_rhs(self):
        initial = initial_state(0.1, 12)
        kappa = historical.LAMBDA / (historical.C_D * historical.L_V)
        gradient = np.array([-kappa * historical.G, -kappa * initial[1], 1 / initial[2]])
        self.assertAlmostEqual(float(gradient @ positive_rhs(0, initial)), 0, places=11)
        result = integrate(initial, t_max=0.01, dt=0.001)
        self.assertLess(float(np.max(np.abs(invariant_residual(result.states)))), 1e-7)


if __name__ == "__main__":
    unittest.main()
