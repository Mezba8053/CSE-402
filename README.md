# PyLBM Vortex-Flowmeter Numerical Project

Newcomers should read [START_HERE.md](START_HERE.md) first.

This all-Python course project uses **PyLBM 0.11.0** as its lattice-Boltzmann solver. Our source code does not implement collision, streaming, D2Q9 populations, MRT matrices or bounce-back formulas manually.

## Solver choice

```text
Library:    PyLBM 0.11.0
Lattice:    D2Q9
Scheme:     Geier central-moment MRT
Backend:    NumPy
Obstacle:   Bouzidi bounce-back
Outlet:     Neumann
Dimension:  2D
```

PyLBM’s official vortex-street example uses the same scheme family. PyLBM performs the LBM time step and returns density and momentum fields.

Our Python code implements the project-specific work:

- triangular/concave/filleted geometry construction;
- case configuration and unit conversion;
- pressure sensor placement;
- vortex-signal extraction;
- peak-to-peak and FFT frequency;
- pressure loss and vorticity;
- trapezoidal and Simpson integration;
- Lagrange, Newton and spline interpolation;
- QR least-squares fitting;
- Golden-Section and Newton optimization;
- result-quality gates and optimum verification.

## Paper relationship

The paper uses C++ OpenLB with a 3D D3Q19 MRT-LES model. This repository uses PyLBM because project code must remain Python. It follows the paper’s problem and workflow but is a reduced 2D educational model, not an exact quantitative reproduction.

## Installation

```bash
python -m pip install -r requirements.txt
```

On Ubuntu/WSL, PyLBM also needs an MPI runtime:

```bash
sudo apt-get update
sudo apt-get install libopenmpi40
```

On older Ubuntu releases the package may be named `libopenmpi3`.

## Commands

Run tests:

```bash
python verify_project.py
```

Run the small case:

```bash
python run_project.py simulate configs/smoke.json
```

Reanalyze saved results:

```bash
python run_project.py analyze data/raw/pylbm_smoke
```

Run a radius sweep and optimization:

```bash
python run_project.py sweep configs/course.json --radii 40 35 30 25 20
```

Verify the fitted optimum with a new flow simulation:

```bash
python run_project.py verify-optimum configs/course.json results/optimization/optimization_summary.json
```

## Data flow

```text
JSON settings
    -> PyLBM geometry and scheme configuration
    -> PyLBM simulation
    -> density and momentum fields
    -> pressure, velocity and vortex signal
    -> frequency, pressure loss and vorticity
    -> geometry-response table
    -> interpolation and QR fit
    -> Golden/Newton optimum
    -> direct verification simulation
```

## Output structure

```text
data/raw/<case>/
|-- config_resolved.json
|-- geometry_summary.json
|-- run_status.json
|-- time_series.csv
|-- metrics.json
|-- field_*.npz
|-- spectrum.npz
|-- vorticity.npy
`-- figures/
    |-- signals.png
    |-- spectrum.png
    `-- flow_fields.png
```

Optimization outputs are stored in `data/processed/` and `results/optimization/`.

## Repository map

| File | Purpose |
|---|---|
| `vortexflow/lbm.py` | Thin PyLBM configuration, execution and sampling adapter |
| `vortexflow/pylbm_compat.py` | Python 3.14 compatibility only; no solver mathematics |
| `vortexflow/geometry.py` | Creates PyLBM triangle/circle geometry elements |
| `vortexflow/postprocess.py` | Signal, pressure and vorticity calculations |
| `vortexflow/numerics.py` | Course numerical algorithms |
| `vortexflow/optimization.py` | Objective, response fitting and searches |
| `vortexflow/pipeline.py` | Runs cases, sweeps and direct verification |
| `vortexflow/plots.py` | Creates figures |
| `configs/` | Reproducible input files |
| `tests/` | PyLBM integration and numerical-method tests |

See [SYSTEM_ARCHITECTURE.md](SYSTEM_ARCHITECTURE.md) for the design and [CODEBASE_GUIDE.md](CODEBASE_GUIDE.md) for the execution path.

## Scientific limits

- The flow domain uses PyLBM’s official vortex-street far-field arrangement, not a detailed no-slip 3D pipe.
- The curved body is constructed by combining a triangle with circular fluid cuts and a circular solid cap.
- A pressure difference across the wake is used as the vortex sensor; momentum-exchange lift is not calculated.
- Course Reynolds numbers are lower than the paper’s methane-pipeline values.
- A full conclusion requires a grid study, a complete radius sweep and direct optimum verification.
