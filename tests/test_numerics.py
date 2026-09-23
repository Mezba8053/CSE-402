from __future__ import annotations

import unittest
import csv
from pathlib import Path
import tempfile

import numpy as np

from vortexflow.numerics import (
    golden_section_maximize,
    lagrange_evaluate,
    natural_cubic_spline,
    newton_coefficients,
    newton_evaluate,
    newton_maximize_polynomial,
    polynomial_evaluate,
    qr_least_squares,
    simpson_uniform,
    trapezoidal,
)
from vortexflow.postprocess import frequency_fft, frequency_peak_to_peak, vorticity
from vortexflow.optimization import optimize_metrics


class NumericalMethodsTests(unittest.TestCase):
    def test_integration(self) -> None:
        x = np.linspace(0.0, 2.0, 101)
        y = x**2
        exact = 8.0 / 3.0
        self.assertAlmostEqual(simpson_uniform(x, y), exact, places=10)
        self.assertLess(abs(trapezoidal(x, y) - exact), 3e-4)

    def test_interpolation(self) -> None:
        x = np.array([-1.0, 0.0, 2.0])
        y = x**2 + 2.0 * x + 3.0
        target = np.linspace(-1.0, 2.0, 20)
        exact = target**2 + 2.0 * target + 3.0
        np.testing.assert_allclose(lagrange_evaluate(x, y, target), exact, atol=1e-12)
        coef = newton_coefficients(x, y)
        np.testing.assert_allclose(newton_evaluate(x, coef, target), exact, atol=1e-12)
        spline = natural_cubic_spline(np.linspace(0, 4, 5), np.linspace(0, 4, 5))
        np.testing.assert_allclose(spline.evaluate(np.linspace(0, 4, 25)), np.linspace(0, 4, 25))

    def test_qr_least_squares(self) -> None:
        x = np.linspace(-2.0, 2.0, 9)
        y = 3.0 - 2.0 * x + 0.5 * x**2
        coefficients, residual = qr_least_squares(x, y, 2)
        np.testing.assert_allclose(coefficients, [3.0, -2.0, 0.5], atol=1e-12)
        self.assertLess(residual, 1e-11)

    def test_optimizers(self) -> None:
        function = lambda value: 5.0 - (value - 2.0) ** 2
        golden = golden_section_maximize(function, -3.0, 6.0)
        self.assertTrue(golden.converged)
        self.assertAlmostEqual(golden.x, 2.0, places=5)
        polynomial = np.array([1.0, 4.0, -1.0])
        newton = newton_maximize_polynomial(polynomial, 0.0, -3.0, 6.0)
        self.assertTrue(newton.converged)
        self.assertAlmostEqual(newton.x, 2.0, places=8)
        self.assertAlmostEqual(float(polynomial_evaluate(polynomial, 2.0)), 5.0)

    def test_signal_frequency(self) -> None:
        time = np.linspace(0.0, 5.0, 1001)
        signal = np.sin(2.0 * np.pi * 7.0 * time)
        peak, _, count = frequency_peak_to_peak(time, signal)
        fft, _, _ = frequency_fft(time, signal)
        self.assertGreater(count, 20)
        self.assertAlmostEqual(peak, 7.0, places=2)
        self.assertLess(abs(fft - 7.0), 1.0 / (time[-1] - time[0]))

    def test_frequency_rejects_roundoff_noise(self) -> None:
        time = np.linspace(0.0, 1.0, 501)
        signal = 1e-18 * np.sin(2.0 * np.pi * 40.0 * time)
        with self.assertRaises(ValueError):
            frequency_peak_to_peak(time, signal)

    def test_vorticity(self) -> None:
        axis = np.linspace(-1.0, 1.0, 31)
        x, y = np.meshgrid(axis, axis)
        omega = vorticity(-y, x, axis[1] - axis[0], axis[1] - axis[0])
        np.testing.assert_allclose(omega[2:-2, 2:-2], 2.0, atol=1e-12)

    def test_optimization_pipeline_writes_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metrics = root / "metrics.csv"
            fields = [
                "case_id", "side_radius_mm", "frequency_hz",
                "lift_amplitude_npm", "mean_pressure_drop_pa_trapezoidal",
                "symmetry_deviation", "frequency_reliable", "stable",
            ]
            with metrics.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for radius, score in [(40, 0.0), (25, 0.8), (35, 0.4), (30, 1.0)]:
                    writer.writerow({
                        "case_id": f"R{radius}", "side_radius_mm": radius,
                        "frequency_hz": 100 + 10 * score,
                        "lift_amplitude_npm": 20 + score,
                        "mean_pressure_drop_pa_trapezoidal": 400 - 40 * score,
                        "symmetry_deviation": 0.10 - 0.02 * score,
                        "frequency_reliable": True, "stable": True,
                    })
            output = root / "optimization"
            summary = optimize_metrics(metrics, output)
            self.assertTrue(summary["golden_section"]["converged"])
            self.assertTrue((output / "response_curves.csv").is_file())
            self.assertTrue((output / "optimization_summary.json").is_file())


if __name__ == "__main__":
    unittest.main()
