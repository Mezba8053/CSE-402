from __future__ import annotations

import csv
import json
import tempfile
import unittest

import numpy as np

from vortexflow.config import SimulationConfig
from vortexflow.geometry import build_pylbm_elements, geometry_summary
from vortexflow.lbm import run_simulation


class PylbmIntegrationTests(unittest.TestCase):
    def test_geometry_is_built_from_library_elements(self) -> None:
        config = SimulationConfig(nx=100, ny=36, steps=60, sample_start=10)
        elements = build_pylbm_elements(config)
        self.assertEqual(len(elements), 4)
        self.assertEqual(elements[0].__class__.__name__, "Triangle")
        self.assertEqual(elements[-1].__class__.__name__, "Circle")
        self.assertGreater(geometry_summary(config)["side_radius_pipe_diameters"], 0)

    def test_short_pylbm_simulation_writes_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            config = SimulationConfig(
                case_id="pylbm_unit_smoke", nx=80, ny=30, steps=50,
                sample_start=20, sample_interval=5, snapshot_interval=50,
                inlet_lattice_velocity=0.03, reynolds_number=60.0,
                output_root=temporary, save_snapshots=True,
            )
            result = run_simulation(config, progress=False)
            self.assertTrue(result.stable)
            self.assertGreater(result.samples, 5)
            self.assertTrue(np.isfinite(result.max_velocity))
            with (config.case_directory / "run_status.json").open(encoding="utf-8") as handle:
                status = json.load(handle)
            self.assertEqual(status["solver"], "PyLBM 0.11.0")
            with (config.case_directory / "time_series.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                fields = csv.DictReader(handle).fieldnames
            self.assertIn("vortex_signal_pa", fields or [])


if __name__ == "__main__":
    unittest.main()
