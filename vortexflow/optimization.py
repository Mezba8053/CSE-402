"""Geometry-response fitting and one-dimensional optimization."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .numerics import (
    golden_section_maximize,
    lagrange_evaluate,
    natural_cubic_spline,
    newton_coefficients,
    newton_evaluate,
    newton_maximize_polynomial,
    polynomial_evaluate,
    qr_least_squares,
)


DEFAULT_WEIGHTS = {
    "frequency": 0.35,
    "lift": 0.20,
    "pressure": 0.30,
    "symmetry": 0.15,
}


def objective_values(rows: list[dict], weights: dict[str, float] | None = None) -> np.ndarray:
    weights = weights or DEFAULT_WEIGHTS
    baseline = rows[0]

    def normalized(row: dict, key: str) -> float:
        denominator = abs(float(baseline[key]))
        if denominator < 1e-14:
            raise ValueError(f"baseline {key} is zero")
        return float(row[key]) / denominator

    return np.asarray(
        [
            weights["frequency"] * normalized(row, "frequency_hz")
            + weights["lift"] * normalized(row, "lift_amplitude_npm")
            - weights["pressure"] * normalized(row, "mean_pressure_drop_pa_trapezoidal")
            - weights["symmetry"] * normalized(row, "symmetry_deviation")
            for row in rows
        ],
        dtype=float,
    )


def optimize_metrics(
    metrics_csv: str | Path,
    output_directory: str | Path,
    weights: dict[str, float] | None = None,
) -> dict:
    with Path(metrics_csv).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < 4:
        raise ValueError("at least four completed geometry cases are required")
    required = [
        "side_radius_mm", "frequency_hz", "lift_amplitude_npm",
        "mean_pressure_drop_pa_trapezoidal", "symmetry_deviation",
    ]
    for row in rows:
        if "frequency_reliable" in row and row["frequency_reliable"].strip().lower() not in {
            "true", "1", "yes"
        }:
            raise ValueError(
                f"unreliable shedding frequency in case {row.get('case_id', '?')}"
            )
        if "stable" in row and row["stable"].strip().lower() not in {"true", "1", "yes"}:
            raise ValueError(f"unstable simulation in case {row.get('case_id', '?')}")
        for name in required:
            value = float(row[name])
            if not np.isfinite(value):
                raise ValueError(f"non-finite {name} in case {row.get('case_id', '?')}")

    # Preserve the first input row as the declared baseline, then sort x and
    # its already-normalized objective together for interpolation.
    objective = objective_values(rows, weights)
    order = np.argsort([float(row["side_radius_mm"]) for row in rows])
    rows = [rows[int(index)] for index in order]
    objective = objective[order]
    x = np.asarray([float(row["side_radius_mm"]) for row in rows])
    if len(np.unique(x)) != len(x):
        raise ValueError("side-radius values must be distinct")
    degree = min(3, len(rows) - 1)
    coefficients, residual = qr_least_squares(x, objective, degree)
    fit = lambda value: float(polynomial_evaluate(coefficients, value))
    golden = golden_section_maximize(fit, float(x[0]), float(x[-1]), tolerance=1e-7)
    newton = newton_maximize_polynomial(
        coefficients, float(np.mean(x)), float(x[0]), float(x[-1])
    )

    newton_coef = newton_coefficients(x, objective)
    spline = natural_cubic_spline(x, objective)
    grid = np.linspace(x[0], x[-1], 401)
    curve = np.column_stack(
        [
            grid,
            lagrange_evaluate(x, objective, grid),
            newton_evaluate(x, newton_coef, grid),
            spline.evaluate(grid),
            polynomial_evaluate(coefficients, grid),
        ]
    )

    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "response_curves.csv",
        curve,
        delimiter=",",
        header="side_radius_mm,lagrange,newton,spline,qr_polynomial",
        comments="",
    )
    summary = {
        "weights": weights or DEFAULT_WEIGHTS,
        "qr_polynomial_coefficients_in_increasing_order": coefficients.tolist(),
        "qr_residual_norm": residual,
        "golden_section": {
            "x_mm": golden.x,
            "objective": golden.value,
            "converged": golden.converged,
            "iterations": len(golden.iterations),
        },
        "newton": {
            "x_mm": newton.x,
            "objective": newton.value,
            "converged": newton.converged,
            "iterations": len(newton.iterations),
        },
        "verification_required": True,
    }
    with (output / "optimization_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    with (output / "golden_history.json").open("w", encoding="utf-8") as handle:
        json.dump(golden.iterations, handle, indent=2)
        handle.write("\n")
    with (output / "newton_history.json").open("w", encoding="utf-8") as handle:
        json.dump(newton.iterations, handle, indent=2)
        handle.write("\n")
    return summary
