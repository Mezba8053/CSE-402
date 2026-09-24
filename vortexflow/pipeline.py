"""Experiment orchestration for simulations, sweeps and verification."""

from __future__ import annotations

from dataclasses import replace
import csv
<<<<<<< HEAD
=======
import json
>>>>>>> 56529ac (added pylbm model for optimized side radius)
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
<<<<<<< HEAD
    "peak_count", "strouhal_number", "lift_amplitude_npm", "symmetry_deviation",
    "mean_pressure_drop_pa_trapezoidal", "mean_pressure_drop_pa_simpson",
    "integration_difference_pa", "max_abs_vorticity_per_s", "stable",
    "mass_relative_change", "max_lattice_velocity",
=======
    "peak_count", "strouhal_number", "signal_amplitude_pa",
    "signal_symmetry_deviation",
    "mean_pressure_drop_pa_trapezoidal", "mean_pressure_drop_pa_simpson",
    "integration_difference_pa", "pressure_stationarity_relative",
    "signal_rms_stationarity_relative", "sampling_stationary",
    "max_abs_vorticity_per_s", "stable",
    "mass_relative_change", "mass_conservation_ok", "max_lattice_velocity",
    "max_lattice_mach", "low_mach_valid", "volumetric_flow_rate_m3ps",
    "mass_flow_rate_kgps",
>>>>>>> 56529ac (added pylbm model for optimized side radius)
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
<<<<<<< HEAD
=======


def verify_optimum(
    base_config_path: str | Path,
    optimization_summary_path: str | Path,
    progress: bool = True,
) -> dict:
    """Run a fresh flow case at the Golden-Section surrogate optimum."""
    summary_path = Path(optimization_summary_path)
    with summary_path.open(encoding="utf-8") as handle:
        summary = json.load(handle)
    radius = float(summary["golden_section"]["x_mm"])
    if not summary["golden_section"].get("converged", False):
        raise ValueError("Golden-Section result did not converge")
    base = load_config(base_config_path)
    tag = f"{radius:.4f}".replace(".", "p")
    config = replace(
        base,
        case_id=f"verified_optimum_R_{tag}mm",
        side_radius_mm=radius,
    )
    metrics = run_case(config, progress=progress)
    target = summary_path.parent / "optimum_verification_metrics.json"
    with target.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2, allow_nan=True)
        handle.write("\n")
    return metrics
>>>>>>> 56529ac (added pylbm model for optimized side radius)
