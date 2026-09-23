"""Command-line interface."""

from __future__ import annotations

import argparse
import json

from .config import load_config
from .optimization import optimize_metrics
from .pipeline import run_case, run_sweep
from .plots import plot_case, plot_optimization
from .postprocess import analyze_case


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
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "simulate":
        metrics = run_case(load_config(args.config), progress=not args.quiet)
        print(json.dumps(metrics, indent=2, allow_nan=True))
    elif args.command == "analyze":
        metrics = analyze_case(args.case_directory)
        plot_case(args.case_directory)
        print(json.dumps(metrics, indent=2, allow_nan=True))
    elif args.command == "sweep":
        metrics, summary = run_sweep(args.config, args.radii, progress=not args.quiet)
        print(f"metrics: {metrics}")
        print(json.dumps(summary, indent=2))
    elif args.command == "optimize":
        summary = optimize_metrics(args.metrics_csv, args.output)
        plot_optimization(args.output, args.metrics_csv)
        print(json.dumps(summary, indent=2))
