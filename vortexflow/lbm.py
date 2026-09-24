"""Thin PyLBM adapter—no collision, streaming, or bounce-back implementation.

The Geier central-moment D2Q9 configuration follows PyLBM's official
``Karman_vortex_street`` example (BSD-3-Clause).  PyLBM owns all LBM numerical
operations.  This module only supplies case parameters, samples returned fields,
and writes the project data contract.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
import time

import numpy as np
import pylbm
import sympy as sp

from .config import SimulationConfig, save_config
from .geometry import build_pylbm_elements, geometry_summary
from .pylbm_compat import apply_python314_boundary_fix


apply_python314_boundary_fix(pylbm)

X, Y = sp.symbols("X Y")
RHO, QX, QY = sp.symbols("rho qx qy")
LA = sp.symbols("lambda", constants=True)


@dataclass(slots=True)
class SimulationResult:
    case_directory: Path
    elapsed_seconds: float
    samples: int
    stable: bool
    max_velocity: float
    mass_relative_change: float


def _inlet_boundary(f, m, x, y, density, velocity):
    """Values requested by PyLBM at the left inlet."""
    m[RHO] = density
    m[QX] = density * velocity
    m[QY] = 0.0


def _geier_scheme(config: SimulationConfig) -> dict:
    """Create the official PyLBM central-moment Navier–Stokes preset."""
    rho0 = 1.0
    velocity = config.inlet_lattice_velocity
    dx = 1.0 / config.ny
    body_height = config.blockage_ratio
    shear_viscosity = rho0 * velocity * body_height / config.reynolds_number
    bulk_viscosity = 10.0 * shear_viscosity
    ux, uy = QX / RHO, QY / RHO
    rho_u2 = RHO * (ux**2 + uy**2)
    polynomials = [
        1, X, Y, X**2 + Y**2, X * Y**2, Y * X**2,
        X**2 * Y**2, X**2 - Y**2, X * Y,
    ]
    equilibrium = [
        RHO, QX, QY,
        rho_u2 + sp.Rational(2, 3) * RHO * LA**2,
        QX * (LA**2 / 3 + uy**2),
        QY * (LA**2 / 3 + ux**2),
        RHO * (LA**2 / 3 + ux**2) * (LA**2 / 3 + uy**2),
        RHO * (ux**2 - uy**2),
        RHO * ux * uy,
    ]
    factor = 3.0 / (rho0 * dx)
    s_bulk = 1.0 / (0.5 + factor * (bulk_viscosity - 2 * shear_viscosity / 3))
    s_shear = 1.0 / (0.5 + factor * shear_viscosity)
    relaxation = [0.0, 0.0, 0.0, s_bulk, s_bulk, s_bulk, s_bulk, s_shear, s_shear]
    width = config.nx / config.ny
    return {
        # Match PyLBM's official vortex-street far-field arrangement: the
        # inlet, top and bottom carry the prescribed stream velocity, while
        # the right boundary is a Neumann outlet.
        "box": {"x": [0.0, width], "y": [0.0, 1.0], "label": [0, 2, 0, 0]},
        "elements": build_pylbm_elements(config),
        "space_step": dx,
        "scheme_velocity": 1.0,
        "schemes": [{
            "velocities": list(range(9)),
            "polynomials": polynomials,
            "relaxation_parameters": relaxation,
            "equilibrium": equilibrium,
            "conserved_moments": [RHO, QX, QY],
        }],
        "parameters": {LA: 1.0},
        "init": {
            RHO: rho0,
            QX: rho0 * velocity,
            QY: rho0 * velocity * config.inlet_perturbation_fraction,
        },
        "boundary_conditions": {
            0: {
                "method": {0: pylbm.bc.BouzidiBounceBack},
                "value": (_inlet_boundary, (rho0, velocity)),
            },
            1: {"method": {0: pylbm.bc.BouzidiBounceBack}},
            2: {"method": {0: pylbm.bc.NeumannX}},
        },
        "generator": "numpy",
        "relative_velocity": [QX / RHO, QY / RHO],
    }


def _fields(solution, config: SimulationConfig) -> tuple[np.ndarray, ...]:
    rho = np.asarray(solution.m[RHO], dtype=float).T
    safe = np.where(np.abs(rho) > 1e-14, rho, 1.0)
    ux_lattice = np.asarray(solution.m[QX], dtype=float).T / safe
    uy_lattice = np.asarray(solution.m[QY], dtype=float).T / safe
    scale = config.physical_velocity_mps / config.inlet_lattice_velocity
    ux = ux_lattice * scale
    uy = uy_lattice * scale
    interior = np.asarray(solution.domain.in_or_out[1:-1, 1:-1]).T
    # PyLBM uses -1 for cells outside the fluid and positive material codes
    # (typically 999) for interior fluid cells.
    solid = interior < 0
    ux[solid] = 0.0
    uy[solid] = 0.0
    return rho, ux_lattice, uy_lattice, ux, uy, solid


def _probe_values(
    rho: np.ndarray, config: SimulationConfig, pressure_scale: float
) -> tuple[float, float, float, float]:
    x_face = int(config.generator_x_ratio * config.nx)
    body_height = max(4, int(config.blockage_ratio * config.ny))
    body_length = max(2, int(config.generator_length_ratio * body_height))
    x_up = max(2, x_face - body_height)
    x_down = min(config.nx - 3, x_face + body_length + 2 * body_height)
    center = config.ny // 2
    band = slice(max(1, center - body_height), min(config.ny - 1, center + body_height))
    p_up = float(np.mean((rho[band, x_up] - 1.0) / 3.0)) * pressure_scale
    p_down = float(np.mean((rho[band, x_down] - 1.0) / 3.0)) * pressure_scale

    sensor_x = min(config.nx - 3, x_face + body_length + body_height)
    offset = max(2, body_height // 2)
    upper = slice(center + offset - 1, center + offset + 2)
    lower = slice(center - offset - 1, center - offset + 2)
    x_slice = slice(sensor_x - 1, sensor_x + 2)
    p_upper = float(np.mean((rho[upper, x_slice] - 1.0) / 3.0)) * pressure_scale
    p_lower = float(np.mean((rho[lower, x_slice] - 1.0) / 3.0)) * pressure_scale
    return p_up, p_down, p_up - p_down, p_upper - p_lower


def run_simulation(config: SimulationConfig, progress: bool = True) -> SimulationResult:
    config.validate()
    output = config.case_directory
    output.mkdir(parents=True, exist_ok=True)
    save_config(config, output / "config_resolved.json")
    with (output / "geometry_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(geometry_summary(config), handle, indent=2)
        handle.write("\n")

    start_clock = time.perf_counter()
    solution = pylbm.Simulation(_geier_scheme(config))
    rho, ux_l, uy_l, ux, uy, solid = _fields(solution, config)
    initial_mass = float(np.sum(rho[~solid]))
    velocity_scale = config.physical_velocity_mps / config.inlet_lattice_velocity
    pressure_scale = config.physical_density_kgpm3 * velocity_scale**2
    max_velocity_seen = 0.0
    samples = 0

    with (output / "time_series.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "step", "time_s", "vortex_signal_pa", "p_upstream_pa",
            "p_downstream_pa", "pressure_drop_pa", "max_lattice_velocity",
        ])
        for step in range(config.steps + 1):
            if step:
                solution.one_time_step()
            needs_fields = (
                step >= config.sample_start and step % config.sample_interval == 0
            ) or (
                config.save_snapshots
                and (step == config.steps or step % config.snapshot_interval == 0)
            ) or step % max(1, config.steps // 10) == 0
            if not needs_fields:
                continue
            rho, ux_l, uy_l, ux, uy, solid = _fields(solution, config)
            speed_lattice = np.sqrt(ux_l**2 + uy_l**2)
            maximum = float(np.max(speed_lattice[~solid]))
            max_velocity_seen = max(max_velocity_seen, maximum)
            if not np.isfinite(rho).all() or not np.isfinite(maximum) or maximum > 0.5:
                raise RuntimeError(
                    f"PyLBM simulation became unstable at step {step}; max velocity={maximum}"
                )
            if step >= config.sample_start and step % config.sample_interval == 0:
                p_up, p_down, pressure_drop, signal = _probe_values(
                    rho, config, pressure_scale
                )
                writer.writerow([
                    step, step * config.physical_dt_s, signal, p_up, p_down,
                    pressure_drop, maximum,
                ])
                samples += 1
            if config.save_snapshots and (
                step == config.steps or step % config.snapshot_interval == 0
            ):
                np.savez_compressed(
                    output / f"field_{step:07d}.npz",
                    rho=rho, ux=ux, uy=uy, solid=solid, body=solid,
                    step=step, time_s=step * config.physical_dt_s,
                )
            if progress and step % max(1, config.steps // 10) == 0:
                print(
                    f"[{config.case_id}] PyLBM step {step:6d}/{config.steps} "
                    f"max|u|={maximum:.5f}"
                )

    final_mass = float(np.sum(rho[~solid]))
    mass_change = abs(final_mass - initial_mass) / initial_mass
    max_mach = float(max_velocity_seen * np.sqrt(3.0))
    status = {
        "solver": "PyLBM 0.11.0",
        "scheme": "D2Q9 Geier central-moment MRT",
        "backend": "NumPy",
        "stable": True,
        "elapsed_seconds": time.perf_counter() - start_clock,
        "samples": samples,
        "max_lattice_velocity": max_velocity_seen,
        "max_lattice_mach": max_mach,
        "low_mach_valid": bool(max_mach <= 0.30),
        "mass_relative_change": mass_change,
        # The domain has open inlet/outlet boundaries, so total domain mass
        # can drift slightly while fluid enters and leaves.  Treat a <=2%
        # drift as acceptable and report the measured value separately.
        "mass_conservation_ok": bool(mass_change <= 0.02),
    }
    with (output / "run_status.json").open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2)
        handle.write("\n")
    return SimulationResult(
        output, status["elapsed_seconds"], samples, True,
        max_velocity_seen, mass_change,
    )
