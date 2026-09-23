"""Configuration loading and validation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path


@dataclass(slots=True)
class SimulationConfig:
    case_id: str = "baseline"
    nx: int = 180
    ny: int = 54
    steps: int = 3500
    sample_start: int = 800
    sample_interval: int = 4
    snapshot_interval: int = 500
    inlet_lattice_velocity: float = 0.08
    inlet_perturbation_fraction: float = 0.01
    reynolds_number: float = 200.0
    smagorinsky_constant: float = 0.12
    blockage_ratio: float = 0.30
    generator_x_ratio: float = 0.32
    generator_length_ratio: float = 0.42
    incoming_angle_deg: float = 180.0
    side_radius_mm: float = 40.0
    fillet_radius_mm: float = 4.0
    pipe_diameter_m: float = 0.09718
    physical_velocity_mps: float = 10.21
    physical_density_kgpm3: float = 27.1721
    characteristic_width_m: float = 0.028
    output_root: str = "data/raw"
    save_snapshots: bool = True

    @property
    def body_height_cells(self) -> float:
        return self.blockage_ratio * (self.ny - 2)

    @property
    def body_length_cells(self) -> float:
        return self.generator_length_ratio * self.body_height_cells

    @property
    def body_x_cells(self) -> float:
        return self.generator_x_ratio * self.nx

    @property
    def physical_dx_m(self) -> float:
        return self.pipe_diameter_m / (self.ny - 2)

    @property
    def physical_dt_s(self) -> float:
        return self.inlet_lattice_velocity * self.physical_dx_m / self.physical_velocity_mps

    @property
    def lattice_viscosity(self) -> float:
        return self.inlet_lattice_velocity * self.body_height_cells / self.reynolds_number

    @property
    def base_relaxation_time(self) -> float:
        return 0.5 + 3.0 * self.lattice_viscosity

    @property
    def case_directory(self) -> Path:
        return Path(self.output_root) / self.case_id

    def validate(self) -> None:
        if not self.case_id or any(token in self.case_id for token in ("..", "/", "\\")):
            raise ValueError("case_id must be a simple directory name")
        if self.nx < 80 or self.ny < 30:
            raise ValueError("grid must be at least 80 x 30")
        if self.steps < 50:
            raise ValueError("steps must be at least 50")
        if not 0 <= self.sample_start < self.steps:
            raise ValueError("sample_start must lie inside the simulation")
        if self.sample_interval < 1 or self.snapshot_interval < 1:
            raise ValueError("output intervals must be positive")
        if not 0.01 <= self.inlet_lattice_velocity <= 0.15:
            raise ValueError("inlet_lattice_velocity should be between 0.01 and 0.15")
        if not 0.0 <= self.inlet_perturbation_fraction <= 0.05:
            raise ValueError("inlet_perturbation_fraction must be between 0 and 0.05")
        if self.reynolds_number <= 0:
            raise ValueError("reynolds_number must be positive")
        if not 0.1 <= self.blockage_ratio <= 0.5:
            raise ValueError("blockage_ratio must be between 0.1 and 0.5")
        if not 0.15 <= self.generator_x_ratio <= 0.65:
            raise ValueError("generator_x_ratio is outside the supported range")
        if self.side_radius_mm <= 0 or self.fillet_radius_mm <= 0:
            raise ValueError("geometry radii must be positive")
        if self.base_relaxation_time <= 0.5005:
            raise ValueError(
                "relaxation time is too close to 0.5; increase resolution, "
                "lattice velocity, or reduce Reynolds number"
            )
        if self.physical_velocity_mps <= 0 or self.pipe_diameter_m <= 0:
            raise ValueError("physical velocity and diameter must be positive")

    def to_dict(self) -> dict:
        result = asdict(self)
        result.update(
            physical_dx_m=self.physical_dx_m,
            physical_dt_s=self.physical_dt_s,
            body_height_cells=self.body_height_cells,
            body_length_cells=self.body_length_cells,
            base_relaxation_time=self.base_relaxation_time,
            model="D2Q9 MRT-LBM with Smagorinsky LES correction",
        )
        return result


def load_config(path: str | Path) -> SimulationConfig:
    with Path(path).open(encoding="utf-8") as handle:
        values = json.load(handle)
    config = SimulationConfig(**values)
    config.validate()
    return config


def save_config(config: SimulationConfig, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as handle:
        json.dump(config.to_dict(), handle, indent=2)
        handle.write("\n")
