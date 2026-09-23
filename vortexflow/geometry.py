"""Parameterized two-dimensional pipe and bluff-body geometry."""

from __future__ import annotations

import numpy as np

from .config import SimulationConfig


def generator_mask(config: SimulationConfig) -> np.ndarray:
    """Return the solid mask for an extruded generator cross-section.

    The sides use a circular-sagitta approximation controlled by the side
    radius. A circular cap represents the downstream fillet. This is an
    explicit engineering approximation because the paper does not publish CAD
    data or a complete coordinate table.
    """

    y, x = np.indices((config.ny, config.nx), dtype=float)
    y_center = 0.5 * (config.ny - 1)
    x_face = config.body_x_cells
    length = config.body_length_cells
    x_tip = x_face + length
    height = config.body_height_cells
    dx_m = config.physical_dx_m
    side_radius = config.side_radius_mm * 1e-3 / dx_m
    fillet = max(0.75, config.fillet_radius_mm * 1e-3 / dx_m)
    fillet = min(fillet, 0.42 * height, 0.42 * length)

    slope = np.tan(np.deg2rad(config.incoming_angle_deg - 180.0))
    local_face = x_face + slope * (y - y_center)
    x_fillet = x_tip - fillet
    usable = np.maximum(x_fillet - local_face, 1e-12)
    s = np.clip((x - local_face) / usable, 0.0, 1.0)
    base = (1.0 - s) * height / 2.0 + s * fillet
    chord = np.sqrt(usable * usable + (height / 2.0 - fillet) ** 2)
    sagitta = np.minimum(chord * chord / (8.0 * max(side_radius, 1e-12)), 0.35 * height)
    curved_half_height = np.maximum(fillet, base - 4.0 * sagitta * s * (1.0 - s))

    cap_dx = x - x_fillet
    cap_half_height = np.sqrt(np.maximum(0.0, fillet * fillet - cap_dx * cap_dx))
    half_height = np.where(x <= x_fillet, curved_half_height, cap_half_height)

    return (
        (x >= local_face)
        & (x <= x_tip)
        & (np.abs(y - y_center) <= half_height)
    )


def build_solid_mask(config: SimulationConfig) -> tuple[np.ndarray, np.ndarray]:
    walls = np.zeros((config.ny, config.nx), dtype=bool)
    walls[0, :] = True
    walls[-1, :] = True
    body = generator_mask(config)
    return walls | body, body
