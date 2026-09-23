"""Experiment orchestration for simulations, sweeps and verification."""

from __future__ import annotations

from dataclasses import replace
import csv
from pathlib import Path

from .config import SimulationConfig, load_config
from .lbm import run_simulation
from .optimization import optimize_metrics
from .plots import plot_case, plot_optimization
from .postprocess import analyze_case


METRIC_COLUMNS = [
    "case_id", "side_radius_mm", "fillet_radius_mm", "incoming_angle_deg",
    "reynolds_number", "frequency_hz", "frequency_std_hz", "fft_frequency_hz",
    "fft_resolution_hz", "frequency_relative_difference", "frequency_reliable",
    "peak_count", "strouhal_number", "lift_amplitude_npm", "symmetry_deviation",
    "mean_pressure_drop_pa_trapezoidal", "mean_pressure_drop_pa_simpson",
    "integration_difference_pa", "max_abs_vorticity_per_s", "stable",
    "mass_relative_change", "max_lattice_velocity",
]


def run_case(config: SimulationConfig, progress: bool = True) -> dict:
    run_simulation(config, progress=progress)
    metrics = analyze_case(config.case_directory)
    plot_case(config.case_directory)
    return metrics


def write_metrics(rows: list[dict], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in METRIC_COLUMNS})
    return target


def run_sweep(
    base_config_path: str | Path,
    radii_mm: list[float],
    progress: bool = True,
) -> tuple[Path, dict]:
    base = load_config(base_config_path)
    rows = []
    for radius in radii_mm:
        tag = str(radius).replace(".", "p")
        config = replace(base, case_id=f"R_{tag}mm", side_radius_mm=radius)
        rows.append(run_case(config, progress=progress))
    metrics_path = write_metrics(rows, Path("data/processed") / "geometry_metrics.csv")
    summary = optimize_metrics(metrics_path, Path("results") / "optimization")
    plot_optimization(Path("results") / "optimization", metrics_path)
    return metrics_path, summary
