from __future__ import annotations

from dataclasses import replace
import tempfile
import unittest

import numpy as np

from vortexflow.config import SimulationConfig
from vortexflow.geometry import build_solid_mask
from vortexflow.lbm import equilibrium, macroscopic, run_simulation


class LbmTests(unittest.TestCase):
    def test_equilibrium_recovers_macroscopic_state(self) -> None:
        rho = np.full((8, 12), 1.02)
        ux = np.full_like(rho, 0.04)
        uy = np.full_like(rho, -0.01)
        distributions = equilibrium(rho, ux, uy)
        recovered = macroscopic(distributions, np.zeros_like(rho, dtype=bool))
        np.testing.assert_allclose(recovered[0], rho, atol=1e-14)
        np.testing.assert_allclose(recovered[1], ux, atol=1e-14)
        np.testing.assert_allclose(recovered[2], uy, atol=1e-14)

    def test_geometry_contains_walls_and_body(self) -> None:
        config = SimulationConfig(nx=100, ny=36, steps=60, sample_start=10)
        solid, body = build_solid_mask(config)
        self.assertTrue(np.all(solid[0]))
        self.assertTrue(np.all(solid[-1]))
        self.assertGreater(np.count_nonzero(body), 10)
        self.assertTrue(np.all(solid[body]))

    def test_short_simulation_remains_finite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = SimulationConfig(
                case_id="unit_smoke", nx=90, ny=32, steps=80,
                sample_start=20, sample_interval=5, snapshot_interval=80,
                inlet_lattice_velocity=0.05, reynolds_number=60.0,
                output_root=temporary,
            )
            result = run_simulation(config, progress=False)
            self.assertTrue(result.stable)
            self.assertGreater(result.samples, 5)
            self.assertTrue(np.isfinite(result.max_velocity))


if __name__ == "__main__":
    unittest.main()
