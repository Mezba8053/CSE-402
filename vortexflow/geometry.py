"""PyLBM geometry objects for the triangular vortex generator."""

from __future__ import annotations

import numpy as np
import pylbm

from .config import SimulationConfig


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
