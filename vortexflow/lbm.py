"""Vectorized D2Q9 MRT lattice-Boltzmann solver with LES correction."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
import time

import numpy as np

from .config import SimulationConfig, save_config
from .geometry import build_solid_mask


C = np.array(
    [[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
     [1, 1], [-1, 1], [-1, -1], [1, -1]],
    dtype=int,
)
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9,
              1 / 36, 1 / 36, 1 / 36, 1 / 36], dtype=float)
OPPOSITE = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=int)

# Standard d'Humieres D2Q9 moment basis.
M = np.array(
    [
        [1, 1, 1, 1, 1, 1, 1, 1, 1],
        [-4, -1, -1, -1, -1, 2, 2, 2, 2],
        [4, -2, -2, -2, -2, 1, 1, 1, 1],
        [0, 1, 0, -1, 0, 1, -1, -1, 1],
        [0, -2, 0, 2, 0, 1, -1, -1, 1],
        [0, 0, 1, 0, -1, 1, 1, -1, -1],
        [0, 0, -2, 0, 2, 1, 1, -1, -1],
        [0, 1, -1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, 1, -1],
    ],
    dtype=float,
)
M_INV = np.linalg.inv(M)


@dataclass(slots=True)
class SimulationResult:
    case_directory: Path
    elapsed_seconds: float
    samples: int
    stable: bool
    max_velocity: float
    mass_relative_change: float


def equilibrium(rho: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    cu = 3.0 * (
        C[:, 0, None, None] * ux[None, :, :]
        + C[:, 1, None, None] * uy[None, :, :]
    )
    u2 = ux * ux + uy * uy
    return W[:, None, None] * rho[None, :, :] * (
        1.0 + cu + 0.5 * cu * cu - 1.5 * u2[None, :, :]
    )


def macroscopic(f: np.ndarray, solid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rho = np.sum(f, axis=0)
    safe_rho = np.where(rho > 1e-14, rho, 1.0)
    ux = np.sum(f * C[:, 0, None, None], axis=0) / safe_rho
    uy = np.sum(f * C[:, 1, None, None], axis=0) / safe_rho
    ux[solid] = 0.0
    uy[solid] = 0.0
    return rho, ux, uy


def _effective_tau(
    f: np.ndarray,
    feq: np.ndarray,
    rho: np.ndarray,
    tau0: float,
    smagorinsky: float,
) -> np.ndarray:
    noneq = f - feq
    qxx = np.sum(noneq * (C[:, 0] ** 2)[:, None, None], axis=0)
    qyy = np.sum(noneq * (C[:, 1] ** 2)[:, None, None], axis=0)
    qxy = np.sum(noneq * (C[:, 0] * C[:, 1])[:, None, None], axis=0)
    norm_q = np.sqrt(qxx * qxx + qyy * qyy + 2.0 * qxy * qxy)
    correction = 18.0 * smagorinsky * smagorinsky * norm_q / np.maximum(rho, 1e-12)
    return 0.5 * (tau0 + np.sqrt(tau0 * tau0 + correction))


def _collide_mrt(
    f: np.ndarray,
    rho: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    config: SimulationConfig,
    fluid: np.ndarray,
) -> np.ndarray:
    feq = equilibrium(rho, ux, uy)
    tau_eff = _effective_tau(
        f, feq, rho, config.base_relaxation_time, config.smagorinsky_constant
    )
    shear_rate = np.clip(1.0 / tau_eff, 0.02, 1.98)

    moments = np.einsum("ab,byx->ayx", M, f, optimize=True)
    moments_eq = np.einsum("ab,byx->ayx", M, feq, optimize=True)
    relaxation = np.empty_like(moments)
    relaxation[0] = 0.0
    relaxation[1] = 1.10
    relaxation[2] = 1.00
    relaxation[3] = 0.0
    relaxation[4] = 1.20
    relaxation[5] = 0.0
    relaxation[6] = 1.20
    relaxation[7] = shear_rate
    relaxation[8] = shear_rate
    post_moments = moments - relaxation * (moments - moments_eq)
    post = np.einsum("ab,byx->ayx", M_INV, post_moments, optimize=True)
    post[:, ~fluid] = f[:, ~fluid]
    return post


def _momentum_exchange(post: np.ndarray, fluid: np.ndarray, solid: np.ndarray) -> tuple[float, float]:
    force_x = 0.0
    force_y = 0.0
    for i, (cx, cy) in enumerate(C):
        if i == 0:
            continue
        destination_solid = np.roll(solid, shift=(-cy, -cx), axis=(0, 1))
        links = fluid & destination_solid
        population = post[i, links]
        force_x += float(np.sum(2.0 * population * cx))
        force_y += float(np.sum(2.0 * population * cy))
    return force_x, force_y


def _stream_with_bounce_back(post: np.ndarray, fluid: np.ndarray, solid: np.ndarray) -> np.ndarray:
    streamed = np.zeros_like(post)
    for i, (cx, cy) in enumerate(C):
        destination_solid = np.roll(solid, shift=(-cy, -cx), axis=(0, 1))
        moving = post[i].copy()
        moving[destination_solid] = 0.0
        streamed[i] += np.roll(moving, shift=(cy, cx), axis=(0, 1))
        links = fluid & destination_solid
        streamed[OPPOSITE[i], links] += post[i, links]
    return streamed


def _apply_boundaries(
    f: np.ndarray,
    target_ux: float,
    target_uy: float,
    solid: np.ndarray,
) -> None:
    # Zou-He velocity inlet at x=0.
    inlet = ~solid[:, 0]
    f0, f2, f4 = f[0, inlet, 0], f[2, inlet, 0], f[4, inlet, 0]
    f3, f6, f7 = f[3, inlet, 0], f[6, inlet, 0], f[7, inlet, 0]
    rho = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - target_ux)
    f[1, inlet, 0] = f3 + (2.0 / 3.0) * rho * target_ux
    f[5, inlet, 0] = (
        f7 + 0.5 * (f4 - f2) + (1.0 / 6.0) * rho * target_ux
        + 0.5 * rho * target_uy
    )
    f[8, inlet, 0] = (
        f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * target_ux
        - 0.5 * rho * target_uy
    )

    # Zou-He constant-density outlet at x=-1.
    outlet = ~solid[:, -1]
    rho_out = 1.0
    ux = -1.0 + (
        f[0, outlet, -1] + f[2, outlet, -1] + f[4, outlet, -1]
        + 2.0 * (f[1, outlet, -1] + f[5, outlet, -1] + f[8, outlet, -1])
    ) / rho_out
    f[3, outlet, -1] = f[1, outlet, -1] - (2.0 / 3.0) * rho_out * ux
    f[7, outlet, -1] = (
        f[5, outlet, -1] + 0.5 * (f[2, outlet, -1] - f[4, outlet, -1])
        - (1.0 / 6.0) * rho_out * ux
    )
    f[6, outlet, -1] = (
        f[8, outlet, -1] + 0.5 * (f[4, outlet, -1] - f[2, outlet, -1])
        - (1.0 / 6.0) * rho_out * ux
    )


def _physical_scales(config: SimulationConfig) -> tuple[float, float]:
    velocity_scale = config.physical_dx_m / config.physical_dt_s
    pressure_scale = config.physical_density_kgpm3 * velocity_scale * velocity_scale
    force_per_depth_scale = pressure_scale * config.physical_dx_m
    return pressure_scale, force_per_depth_scale


def run_simulation(config: SimulationConfig, progress: bool = True) -> SimulationResult:
    config.validate()
    output = config.case_directory
    output.mkdir(parents=True, exist_ok=True)
    save_config(config, output / "config_resolved.json")

    solid, body = build_solid_mask(config)
    fluid = ~solid
    rho = np.ones((config.ny, config.nx), dtype=float)
    ux = np.zeros_like(rho)
    uy = np.zeros_like(rho)
    f = equilibrium(rho, ux, uy)
    initial_mass = float(np.sum(rho[fluid]))
    pressure_scale, force_scale = _physical_scales(config)

    time_series_path = output / "time_series.csv"
    start_clock = time.perf_counter()
    sample_count = 0
    stable = True
    max_velocity_seen = 0.0

    with time_series_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["step", "time_s", "lift_npm", "drag_npm", "p_upstream_pa",
             "p_downstream_pa", "pressure_drop_pa", "max_lattice_velocity"]
        )

        for step in range(config.steps + 1):
            rho, ux, uy = macroscopic(f, solid)
            speed = np.sqrt(ux * ux + uy * uy)
            max_velocity = float(np.max(speed[fluid]))
            max_velocity_seen = max(max_velocity_seen, max_velocity)
            if not np.isfinite(f).all() or not np.isfinite(max_velocity) or max_velocity > 0.5:
                stable = False
                raise RuntimeError(
                    f"simulation became unstable at step {step}; max lattice velocity={max_velocity}"
                )

            post = _collide_mrt(f, rho, ux, uy, config, fluid)
            drag_lattice, lift_lattice = _momentum_exchange(post, fluid, body)
            f = _stream_with_bounce_back(post, fluid, solid)

            ramp_steps = max(1, int(0.15 * config.steps))
            ramp = min(1.0, step / ramp_steps)
            ramp = ramp * ramp * (3.0 - 2.0 * ramp)
            # A perfectly symmetric discrete problem can remain on the unstable
            # symmetric solution forever.  A smooth, startup-only transverse
            # pulse seeds the physical antisymmetric shedding mode without
            # forcing the mature wake at an imposed frequency.
            perturbation = 0.0
            if step <= ramp_steps:
                perturbation = (
                    config.inlet_perturbation_fraction
                    * config.inlet_lattice_velocity
                    * np.sin(np.pi * step / ramp_steps)
                )
            _apply_boundaries(
                f,
                ramp * config.inlet_lattice_velocity,
                perturbation,
                solid,
            )

            if step >= config.sample_start and step % config.sample_interval == 0:
                rho_now, _, _ = macroscopic(f, solid)
                x_up = max(2, int(config.body_x_cells - 1.5 * config.body_height_cells))
                x_down = min(config.nx - 3, int(
                    config.body_x_cells + config.body_length_cells
                    + 1.5 * config.body_height_cells
                ))
                center_band = slice(config.ny // 3, 2 * config.ny // 3)
                p_up = float(np.mean((rho_now[center_band, x_up] - 1.0) / 3.0)) * pressure_scale
                p_down = float(np.mean((rho_now[center_band, x_down] - 1.0) / 3.0)) * pressure_scale
                writer.writerow(
                    [step, step * config.physical_dt_s,
                     lift_lattice * force_scale, drag_lattice * force_scale,
                     p_up, p_down, p_up - p_down, max_velocity]
                )
                sample_count += 1

            if config.save_snapshots and (
                step == config.steps or step % config.snapshot_interval == 0
            ):
                np.savez_compressed(
                    output / f"field_{step:07d}.npz",
                    rho=rho,
                    ux=ux * config.physical_velocity_mps / config.inlet_lattice_velocity,
                    uy=uy * config.physical_velocity_mps / config.inlet_lattice_velocity,
                    solid=solid,
                    body=body,
                    step=step,
                    time_s=step * config.physical_dt_s,
                )

            if progress and step % max(1, config.steps // 10) == 0:
                print(
                    f"[{config.case_id}] step {step:6d}/{config.steps} "
                    f"max|u|={max_velocity:.5f}"
                )

    rho_final, _, _ = macroscopic(f, solid)
    final_mass = float(np.sum(rho_final[fluid]))
    mass_change = abs(final_mass - initial_mass) / initial_mass
    elapsed = time.perf_counter() - start_clock
    status = {
        "stable": stable,
        "elapsed_seconds": elapsed,
        "samples": sample_count,
        "max_lattice_velocity": max_velocity_seen,
        "mass_relative_change": mass_change,
    }
    with (output / "run_status.json").open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2)
        handle.write("\n")
    return SimulationResult(output, elapsed, sample_count, stable, max_velocity_seen, mass_change)
