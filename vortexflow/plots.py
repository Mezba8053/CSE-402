"""Reproducible plots generated from saved result files."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _save_figure(fig: plt.Figure, target: Path) -> None:
    """Save repeatably on Windows/WSL, where in-place PNG overwrite can fail."""
    target.unlink(missing_ok=True)
    fig.savefig(target, dpi=160)


def plot_case(case_directory: str | Path) -> list[Path]:
    case = Path(case_directory)
    figures = case / "figures"
    figures.mkdir(exist_ok=True)
    with (case / "time_series.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    time = np.asarray([float(row["time_s"]) for row in rows])
<<<<<<< HEAD
    lift = np.asarray([float(row["lift_npm"]) for row in rows])
=======
    signal = np.asarray([float(row["vortex_signal_pa"]) for row in rows])
>>>>>>> 56529ac (added pylbm model for optimized side radius)
    pressure = np.asarray([float(row["pressure_drop_pa"]) for row in rows])

    created: list[Path] = []
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
<<<<<<< HEAD
    axes[0].plot(time, lift, linewidth=1.0)
    axes[0].set_ylabel("Lift per depth (N/m)")
=======
    axes[0].plot(time, signal, linewidth=1.0)
    axes[0].set_ylabel("Vortex sensor Δp (Pa)")
>>>>>>> 56529ac (added pylbm model for optimized side radius)
    axes[0].grid(alpha=0.3)
    axes[1].plot(time, pressure, linewidth=1.0, color="tab:red")
    axes[1].set_xlabel("Physical time (s)")
    axes[1].set_ylabel("Pressure drop (Pa)")
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    target = figures / "signals.png"
    _save_figure(fig, target)
    plt.close(fig)
    created.append(target)

    spectrum = np.load(case / "spectrum.npz")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(spectrum["frequency"], spectrum["amplitude"], linewidth=1.0)
<<<<<<< HEAD
=======
    if len(spectrum["frequency"]) > 1:
        dominant_index = int(np.argmax(spectrum["amplitude"][1:]) + 1)
        dominant_frequency = float(spectrum["frequency"][dominant_index])
        upper = min(float(spectrum["frequency"][-1]), 5.0 * dominant_frequency)
        if upper > 0:
            ax.set_xlim(0.0, upper)
>>>>>>> 56529ac (added pylbm model for optimized side radius)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("FFT amplitude")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    target = figures / "spectrum.png"
    _save_figure(fig, target)
    plt.close(fig)
    created.append(target)

    snapshots = sorted(case.glob("field_*.npz"))
    if snapshots:
        field = np.load(snapshots[-1])
        speed = np.sqrt(field["ux"] ** 2 + field["uy"] ** 2)
        speed = np.ma.masked_where(field["solid"], speed)
        omega = np.load(case / "vorticity.npy")
        fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
        image = axes[0].imshow(speed, origin="lower", aspect="auto", cmap="viridis")
        axes[0].set_ylabel("y cell")
        axes[0].set_title("Velocity magnitude (m/s)")
        fig.colorbar(image, ax=axes[0])
        limit = np.nanpercentile(np.abs(omega), 98)
        vort = axes[1].imshow(
            omega, origin="lower", aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit
        )
        axes[1].set_xlabel("x cell")
        axes[1].set_ylabel("y cell")
        axes[1].set_title("Vorticity (1/s)")
        fig.colorbar(vort, ax=axes[1])
        fig.tight_layout()
        target = figures / "flow_fields.png"
        _save_figure(fig, target)
        plt.close(fig)
        created.append(target)
    return created


def plot_optimization(output_directory: str | Path, metrics_csv: str | Path) -> Path:
    output = Path(output_directory)
    curve = np.genfromtxt(output / "response_curves.csv", delimiter=",", names=True)
    with Path(metrics_csv).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (output / "optimization_summary.json").open(encoding="utf-8") as handle:
        summary = json.load(handle)
    x = np.asarray([float(row["side_radius_mm"]) for row in rows])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(curve["side_radius_mm"], curve["lagrange"], label="Lagrange", alpha=0.75)
    ax.plot(curve["side_radius_mm"], curve["newton"], "--", label="Newton interpolation")
    ax.plot(curve["side_radius_mm"], curve["spline"], label="Natural spline")
    ax.plot(curve["side_radius_mm"], curve["qr_polynomial"], linewidth=2, label="QR fit")
    ax.axvline(summary["golden_section"]["x_mm"], color="black", linestyle=":", label="Golden optimum")
    ax.set_xlabel("Side radius (mm)")
    ax.set_ylabel("Normalized objective")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    target = output / "optimization.png"
    _save_figure(fig, target)
    plt.close(fig)
    return target
