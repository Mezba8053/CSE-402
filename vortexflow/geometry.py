<<<<<<< HEAD
"""Parameterized two-dimensional pipe and bluff-body geometry."""
=======
"""PyLBM geometry objects for the triangular vortex generator."""
>>>>>>> 56529ac (added pylbm model for optimized side radius)

from __future__ import annotations

import numpy as np
<<<<<<< HEAD
=======
import pylbm
>>>>>>> 56529ac (added pylbm model for optimized side radius)

from .config import SimulationConfig


<<<<<<< HEAD
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
=======
def _arc_center(start: np.ndarray, end: np.ndarray, radius: float, upper: bool) -> np.ndarray:
    chord_vector = end - start
    chord = float(np.linalg.norm(chord_vector))
    if radius <= 0.5 * chord:
        radius = 0.5 * chord + 1e-6
    midpoint = 0.5 * (start + end)
    normal = np.array([-chord_vector[1], chord_vector[0]]) / chord
    if (normal[1] > 0) != upper:
        normal *= -1.0
    distance = np.sqrt(max(radius * radius - 0.25 * chord * chord, 0.0))
    return midpoint + distance * normal


def build_pylbm_elements(config: SimulationConfig) -> list:
    """Return library geometry elements controlled by the paper parameters.

    A solid triangle supplies the basic bluff body.  Two fluid circles cut
    concave side arcs and one solid circle supplies the rounded downstream cap.
    PyLBM performs the actual grid/domain construction and wall distances.
    """

    width = config.nx / config.ny
    center_y = 0.5
    x_face = config.generator_x_ratio * width
    height = config.blockage_ratio
    length = config.generator_length_ratio * height
    x_tip = x_face + length
    side_radius = config.side_radius_mm * 1e-3 / config.pipe_diameter_m
    fillet = config.fillet_radius_mm * 1e-3 / config.pipe_diameter_m
    fillet = min(fillet, 0.40 * height, 0.40 * length)

    top = np.array([x_face, center_y + 0.5 * height])
    bottom = np.array([x_face, center_y - 0.5 * height])
    cap_top = np.array([x_tip - fillet, center_y + fillet])
    cap_bottom = np.array([x_tip - fillet, center_y - fillet])
    upper_center = _arc_center(top, cap_top, side_radius, upper=True)
    lower_center = _arc_center(bottom, cap_bottom, side_radius, upper=False)

    return [
        pylbm.Triangle(
            top,
            bottom - top,
            np.array([x_tip, center_y]) - top,
            label=1,
            isfluid=False,
        ),
        pylbm.Circle(upper_center, side_radius, label=1, isfluid=True),
        pylbm.Circle(lower_center, side_radius, label=1, isfluid=True),
        pylbm.Circle([x_tip - fillet, center_y], fillet, label=1, isfluid=False),
    ]


def geometry_summary(config: SimulationConfig) -> dict[str, float]:
    return {
        "domain_width_pipe_diameters": config.nx / config.ny,
        "generator_height_pipe_diameters": config.blockage_ratio,
        "generator_length_pipe_diameters": (
            config.generator_length_ratio * config.blockage_ratio
        ),
        "side_radius_pipe_diameters": (
            config.side_radius_mm * 1e-3 / config.pipe_diameter_m
        ),
        "fillet_radius_pipe_diameters": (
            config.fillet_radius_mm * 1e-3 / config.pipe_diameter_m
        ),
    }
>>>>>>> 56529ac (added pylbm model for optimized side radius)
