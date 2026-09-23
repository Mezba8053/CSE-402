"""Signal, flow-field and case-level post-processing."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from .numerics import simpson_uniform, trapezoidal


def _linear_detrend(time: np.ndarray, signal: np.ndarray) -> np.ndarray:
    """Remove offset and linear drift without requiring SciPy."""
    centered_time = time - np.mean(time)
    design = np.column_stack((np.ones_like(centered_time), centered_time))
    coefficients, *_ = np.linalg.lstsq(design, signal, rcond=None)
    return signal - design @ coefficients


def detect_peaks(
    time: np.ndarray,
    signal: np.ndarray,
    minimum_distance: int = 3,
    minimum_height: float | None = None,
) -> np.ndarray:
    if len(signal) < 3:
        return np.array([], dtype=int)
    threshold = float(np.mean(signal)) if minimum_height is None else minimum_height
    candidates = np.flatnonzero(
        (signal[1:-1] > signal[:-2])
        & (signal[1:-1] >= signal[2:])
        & (signal[1:-1] > threshold)
    ) + 1
    selected: list[int] = []
    for index in candidates:
        if not selected or index - selected[-1] >= minimum_distance:
            selected.append(int(index))
        elif signal[index] > signal[selected[-1]]:
            selected[-1] = int(index)
    return np.asarray(selected, dtype=int)


def frequency_peak_to_peak(time: np.ndarray, signal: np.ndarray) -> tuple[float, float, int]:
    if len(time) != len(signal) or len(time) < 8:
        raise ValueError("time and signal must have at least eight matching samples")
    detrended = _linear_detrend(time, signal)
    fft_guess, _, _ = frequency_fft(time, signal)
    sample_period = float(np.mean(np.diff(time)))
    samples_per_cycle = 1.0 / max(fft_guess * sample_period, 1e-15)
    minimum_distance = max(3, int(0.55 * samples_per_cycle))
    amplitude = float(np.sqrt(np.mean(detrended * detrended)))
    if amplitude <= 100.0 * np.finfo(float).eps * max(1.0, float(np.max(np.abs(signal)))):
        raise ValueError("signal amplitude is indistinguishable from round-off noise")
    peaks = detect_peaks(
        time,
        detrended,
        minimum_distance=minimum_distance,
        minimum_height=0.15 * amplitude,
    )
    if len(peaks) < 3:
        raise ValueError("at least three peaks are required for a reliable period estimate")
    periods = np.diff(time[peaks])
    median = float(np.median(periods))
    accepted = periods[np.abs(periods - median) <= 0.25 * median]
    if len(accepted) < 2:
        raise ValueError("peak periods are not sufficiently consistent")
    frequencies = 1.0 / accepted
    estimate = float(1.0 / np.mean(accepted))
    uncertainty = float(np.std(frequencies, ddof=1))
    if uncertainty / estimate > 0.20:
        raise ValueError("peak-to-peak frequency has more than 20% relative scatter")
    return estimate, uncertainty, len(peaks)


def frequency_fft(time: np.ndarray, signal: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    spacing = np.diff(time)
    if not np.allclose(spacing, spacing[0], rtol=1e-5, atol=1e-12):
        raise ValueError("FFT requires uniformly sampled data")
    centered = _linear_detrend(time, signal)
    window = np.hanning(len(centered))
    spectrum = np.abs(np.fft.rfft(centered * window))
    frequencies = np.fft.rfftfreq(len(centered), d=spacing[0])
    if len(spectrum) < 2:
        raise ValueError("not enough samples for FFT")
    # Ignore frequencies that cannot complete at least two cycles in the
    # observation window; those bins mostly represent residual startup drift.
    duration = float(time[-1] - time[0])
    first_valid = max(1, int(np.searchsorted(frequencies, 2.0 / duration)))
    if first_valid >= len(spectrum):
        raise ValueError("time window is too short for frequency estimation")
    index = int(np.argmax(spectrum[first_valid:]) + first_valid)
    return float(frequencies[index]), frequencies, spectrum


def vorticity(ux: np.ndarray, uy: np.ndarray, dx: float, dy: float) -> np.ndarray:
    if ux.shape != uy.shape or ux.ndim != 2:
        raise ValueError("ux and uy must be two-dimensional arrays with matching shapes")
    dv_dx = np.gradient(uy, dx, axis=1, edge_order=2)
    du_dy = np.gradient(ux, dy, axis=0, edge_order=2)
    return dv_dx - du_dy


def analyze_case(case_directory: str | Path) -> dict:
    case = Path(case_directory)
    with (case / "config_resolved.json").open(encoding="utf-8") as handle:
        config = json.load(handle)
    with (case / "run_status.json").open(encoding="utf-8") as handle:
        status = json.load(handle)

    columns: dict[str, list[float]] = {}
    with (case / "time_series.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for name in reader.fieldnames or []:
            columns[name] = []
        for row in reader:
            for name in columns:
                columns[name].append(float(row[name]))
    arrays = {name: np.asarray(values, dtype=float) for name, values in columns.items()}
    if len(arrays.get("time_s", [])) < 4:
        raise ValueError("case does not contain enough time-series samples")

    time = arrays["time_s"]
    lift = arrays["lift_npm"]
    pressure_drop = arrays["pressure_drop_pa"]
    frequency = float("nan")
    frequency_std = float("nan")
    peak_count = 0
    try:
        frequency, frequency_std, peak_count = frequency_peak_to_peak(time, lift)
    except ValueError:
        pass
    fft_frequency, fft_axis, fft_amplitude = frequency_fft(time, lift)
    fft_resolution = 1.0 / (time[-1] - time[0])
    relative_difference = (
        abs(frequency - fft_frequency) / fft_frequency
        if np.isfinite(frequency) and fft_frequency > 0 else float("nan")
    )
    frequency_reliable = bool(
        np.isfinite(frequency)
        and abs(frequency - fft_frequency) <= max(2.0 * fft_resolution, 0.20 * fft_frequency)
    )
    duration = time[-1] - time[0]
    pressure_trapezoidal = trapezoidal(time, pressure_drop) / duration
    pressure_simpson = float("nan")
    if (len(time) - 1) % 2 == 0:
        pressure_simpson = simpson_uniform(time, pressure_drop) / duration

    lift_max = float(np.max(lift))
    lift_min = float(np.min(lift))
    denominator = abs(lift_max - lift_min)
    symmetry = abs(lift_max + lift_min) / denominator if denominator else float("nan")
    strouhal = (
        frequency * config["characteristic_width_m"] / config["physical_velocity_mps"]
        if np.isfinite(frequency) else float("nan")
    )

    snapshots = sorted(case.glob("field_*.npz"))
    if snapshots:
        field = np.load(snapshots[-1])
        omega = vorticity(
            field["ux"], field["uy"], config["physical_dx_m"], config["physical_dx_m"]
        )
        omega = np.where(field["solid"], np.nan, omega)
        np.save(case / "vorticity.npy", omega)
        max_abs_vorticity = float(np.nanmax(np.abs(omega)))
    else:
        max_abs_vorticity = float("nan")

    metrics = {
        "case_id": config["case_id"],
        "side_radius_mm": config["side_radius_mm"],
        "fillet_radius_mm": config["fillet_radius_mm"],
        "incoming_angle_deg": config["incoming_angle_deg"],
        "reynolds_number": config["reynolds_number"],
        "frequency_hz": frequency,
        "frequency_std_hz": frequency_std,
        "fft_frequency_hz": fft_frequency,
        "fft_resolution_hz": fft_resolution,
        "frequency_relative_difference": relative_difference,
        "frequency_reliable": frequency_reliable,
        "peak_count": peak_count,
        "strouhal_number": strouhal,
        "lift_amplitude_npm": 0.5 * denominator,
        "symmetry_deviation": symmetry,
        "mean_pressure_drop_pa_trapezoidal": pressure_trapezoidal,
        "mean_pressure_drop_pa_simpson": pressure_simpson,
        "integration_difference_pa": abs(pressure_trapezoidal - pressure_simpson)
        if np.isfinite(pressure_simpson) else float("nan"),
        "max_abs_vorticity_per_s": max_abs_vorticity,
        "stable": status["stable"],
        "mass_relative_change": status["mass_relative_change"],
        "max_lattice_velocity": status["max_lattice_velocity"],
    }
    with (case / "metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2, allow_nan=True)
        handle.write("\n")
    np.savez_compressed(case / "spectrum.npz", frequency=fft_axis, amplitude=fft_amplitude)
    return metrics
