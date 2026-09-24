"""Command-line interface."""

from __future__ import annotations

import argparse
<<<<<<< HEAD
import json

from .config import load_config
from .optimization import optimize_metrics
from .pipeline import run_case, run_sweep
=======
import csv
import json

from .config import load_config
from .optimization import optimize_metrics, objective_values
from .pipeline import run_case, run_sweep, verify_optimum
>>>>>>> 56529ac (added pylbm model for optimized side radius)
from .plots import plot_case, plot_optimization
from .postprocess import analyze_case


<<<<<<< HEAD
=======
def print_numerical_summary(metrics: dict) -> None:
    """Print the numerical-method outputs most useful for a course demo."""
    radius = metrics.get("side_radius_mm", "fixed by configuration")
    print("\nNUMERICAL METHODS RESULT (execution order)")
    print(f"  Case geometry: side radius = {radius} mm")
    print(f"  1. LBM discretization: {metrics.get('solver', 'PyLBM')} (D2Q9 central-moment MRT)")
    print(f"  2. Stability/error checks: stable={metrics.get('stable')}, mass conservation={metrics.get('mass_conservation_ok')}, low-Mach={metrics.get('low_mach_valid')}")
    print(f"     mass drift={100.0 * float(metrics.get('mass_relative_change', float('nan'))):.3f}%")
    print(f"     flow rate: Q={metrics.get('volumetric_flow_rate_m3ps', float('nan')):.6f} m^3/s, m_dot={metrics.get('mass_flow_rate_kgps', float('nan')):.6f} kg/s")
    print(f"  3. Time-series sampling: stationary={metrics.get('sampling_stationary')}")
    print(f"  4. FFT frequency analysis: {metrics.get('fft_frequency_hz', float('nan')):.4f} Hz (resolution {metrics.get('fft_resolution_hz', float('nan')):.4f} Hz)")
    print(f"  5. Peak interpolation/frequency: {metrics.get('frequency_hz', float('nan')):.4f} Hz")
    print(f"  6. Strouhal calculation: {metrics.get('strouhal_number', float('nan')):.6f}")
    print(f"  7. Trapezoidal integration: {metrics.get('mean_pressure_drop_pa_trapezoidal', float('nan')):.6f} Pa")
    print(f"  8. Simpson integration: {metrics.get('mean_pressure_drop_pa_simpson', float('nan')):.6f} Pa")
    print(f"  9. Integration difference/error: {metrics.get('integration_difference_pa', float('nan')):.6g} Pa")
    print("  10. Interpolation/curve fitting: run with `sweep` (requires multiple radii)")
    print("  11. Golden-Section/Newton optimization: run with `sweep`, then `verify-optimum`")


def print_optimization_summary(summary: dict) -> None:
    """Print interpolation and optimization results with their radius."""
    print("\nINTERPOLATION AND OPTIMIZATION RESULT")
    print(f"  Objective weights: {summary.get('weights', {})}")
    golden = summary.get("golden_section", {})
    newton = summary.get("newton", {})
    print(f"  Golden-Section optimum: radius={golden.get('x_mm', float('nan')):.4f} mm, score={golden.get('objective', float('nan')):.6f}, converged={golden.get('converged')}")
    print(f"  Newton optimum: radius={newton.get('x_mm', float('nan')):.4f} mm, score={newton.get('objective', float('nan')):.6f}, converged={newton.get('converged')}")
    print("  Response curves saved: Lagrange, Newton polynomial, natural cubic spline, QR least-squares")
    print("  Note: these are surrogate-curve predictions; verify-optimum runs a fresh LBM case.")


>>>>>>> 56529ac (added pylbm model for optimized side radius)
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Triangular vortex-flowmeter numerical project")
    commands = parser.add_subparsers(dest="command", required=True)

    simulate = commands.add_parser("simulate", help="run one LBM case and analyze it")
    simulate.add_argument("config")
    simulate.add_argument("--quiet", action="store_true")

    analyze = commands.add_parser("analyze", help="reanalyze an existing case")
    analyze.add_argument("case_directory")

    sweep = commands.add_parser("sweep", help="run a side-radius geometry sweep")
    sweep.add_argument("config")
    sweep.add_argument("--radii", nargs="+", type=float, required=True)
    sweep.add_argument("--quiet", action="store_true")

    optimize = commands.add_parser("optimize", help="optimize an existing metrics table")
    optimize.add_argument("metrics_csv")
    optimize.add_argument("--output", default="results/optimization")
<<<<<<< HEAD
=======

    verify = commands.add_parser(
        "verify-optimum", help="directly simulate a predicted optimum"
    )
    verify.add_argument("config")
    verify.add_argument("optimization_summary")
    verify.add_argument("--quiet", action="store_true")
>>>>>>> 56529ac (added pylbm model for optimized side radius)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "simulate":
        metrics = run_case(load_config(args.config), progress=not args.quiet)
<<<<<<< HEAD
=======
        print_numerical_summary(metrics)
>>>>>>> 56529ac (added pylbm model for optimized side radius)
        print(json.dumps(metrics, indent=2, allow_nan=True))
    elif args.command == "analyze":
        metrics = analyze_case(args.case_directory)
        plot_case(args.case_directory)
<<<<<<< HEAD
        print(json.dumps(metrics, indent=2, allow_nan=True))
    elif args.command == "sweep":
        metrics, summary = run_sweep(args.config, args.radii, progress=not args.quiet)
        print(f"metrics: {metrics}")
=======
        print_numerical_summary(metrics)
        print(json.dumps(metrics, indent=2, allow_nan=True))
    elif args.command == "sweep":
        metrics, summary = run_sweep(args.config, args.radii, progress=not args.quiet)
        print("\nRADIUS SWEEP RESULTS")
        print(f"  Metrics table: {metrics}")
        with metrics.open(newline="", encoding="utf-8") as handle:
            sweep_rows = list(csv.DictReader(handle))
        scores = objective_values(sweep_rows)
        base_config = load_config(args.config)
        nominal_flow_rate = base_config.physical_velocity_mps * 0.25 * 3.141592653589793 * base_config.pipe_diameter_m**2
        for row, score in zip(sweep_rows, scores):
            flow_text = row.get("volumetric_flow_rate_m3ps", "")
            flow_rate = float(flow_text) if flow_text else nominal_flow_rate
            print(
                f"  radius={float(row['side_radius_mm']):.3f} mm | "
                f"f={float(row['frequency_hz']):.4f} Hz | "
                f"FFT={float(row['fft_frequency_hz']):.4f} Hz | "
                f"trap={float(row['mean_pressure_drop_pa_trapezoidal']):.3f} Pa | "
                f"Simpson={float(row['mean_pressure_drop_pa_simpson']):.3f} Pa | "
                f"int.err={float(row['integration_difference_pa']):.3g} Pa | "
                f"St={float(row['strouhal_number']):.6f} | "
                f"J={score:.6f} | mass drift={100*float(row['mass_relative_change']):.2f}% | "
                f"Q={flow_rate:.5f} m^3/s | "
                f"stationary={row['sampling_stationary']}"
            )
        print_optimization_summary(summary)
>>>>>>> 56529ac (added pylbm model for optimized side radius)
        print(json.dumps(summary, indent=2))
    elif args.command == "optimize":
        summary = optimize_metrics(args.metrics_csv, args.output)
        plot_optimization(args.output, args.metrics_csv)
<<<<<<< HEAD
        print(json.dumps(summary, indent=2))
=======
        print_optimization_summary(summary)
        print(json.dumps(summary, indent=2))
    elif args.command == "verify-optimum":
        metrics = verify_optimum(
            args.config, args.optimization_summary, progress=not args.quiet
        )
        print_numerical_summary(metrics)
        print(json.dumps(metrics, indent=2, allow_nan=True))
>>>>>>> 56529ac (added pylbm model for optimized side radius)
