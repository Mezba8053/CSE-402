# Start Here — The Project in Simple Language

## 1. The project in one sentence

We use the Python package **PyLBM** to simulate fluid passing a triangular object, measure the repeating vortex signal, and use Python numerical methods to find a better object shape.

## 2. What is a vortex flowmeter?

A solid object is placed in moving fluid. Alternating swirls form behind it:

```text
flow --->    triangle       clockwise / anticlockwise vortices
                 >              ↻   ↺   ↻   ↺
```

These vortices create a repeating pressure signal. The repetition rate is the vortex frequency. A vortex flowmeter relates this frequency to flow velocity.

We want:

- a clear, regular vortex signal;
- a high enough frequency;
- low pressure loss;
- a geometrically practical vortex generator.

## 3. What PyLBM does for us

We do **not** program the lattice-Boltzmann solver ourselves.

PyLBM supplies:

- the D2Q9 lattice;
- the Geier central-moment MRT scheme;
- distribution functions;
- equilibrium calculations;
- collision;
- streaming;
- Bouzidi bounce-back at the obstacle;
- inlet/far-field boundary treatment;
- Neumann outlet treatment;
- the simulation time-step operation;
- density and momentum fields;
- NumPy execution.

Our code never manually loops over nine lattice directions and never contains a hand-written MRT matrix.

## 4. What we still implement

Using a library does not complete the project automatically. Our Python code supplies:

- the triangular, concave and rounded generator geometry;
- all physical and numerical settings;
- the locations of pressure sensors;
- conversion from lattice units to physical units;
- result files and plots;
- vortex-frequency calculation;
- FFT confirmation;
- pressure-loss calculation;
- vorticity calculation;
- trapezoidal and Simpson integration;
- Lagrange, Newton and spline interpolation;
- QR least-squares fitting;
- Golden-Section Search;
- Newton optimization;
- quality checks and direct optimum verification.

This is where the course’s numerical-methods work is concentrated.

## 5. Why PyLBM instead of OpenLB?

The research paper uses OpenLB, D3Q19, MRT and LES in three dimensions. OpenLB is a C++ framework.

This project must keep its own code in Python, so it uses PyLBM. PyLBM is controlled from Python and provides a documented D2Q9 Von Kármán vortex-street example, which closely matches flow around a bluff body.

The selected library method is:

```text
PyLBM 0.11.0
D2Q9
Geier central-moment multiple-relaxation-time scheme
NumPy backend
Bouzidi obstacle bounce-back
Neumann outlet
```

This is a two-dimensional educational approximation, not an exact reproduction of the paper’s three-dimensional methane-pipe simulation.

## 6. The entire workflow

```text
configs/smoke.json
        |
        v
Create PyLBM triangle + curved cuts + rounded tip
        |
        v
Give geometry and flow settings to PyLBM
        |
        v
PyLBM performs all LBM calculations
        |
        v
PyLBM returns density and momentum fields
        |
        +--> pressure before and after the object
        |
        +--> pressure difference above and below the wake
        |
        `--> velocity field
                 |
                 v
        Our Python post-processing
                 |
      +----------+-----------+
      |          |           |
      v          v           v
 frequency   pressure     vorticity
 peak + FFT  integration  derivatives
      |          |           |
      +----------+-----------+
                 |
                 v
       repeat for several radii
                 |
                 v
 interpolation + QR curve fitting
                 |
                 v
 Golden-Section + Newton optimization
                 |
                 v
 simulate the predicted optimum again
```

## 7. Files to read, in order

1. `configs/smoke.json` — simulation settings.
2. `run_project.py` — starts the program.
3. `vortexflow/pipeline.py` — connects all stages.
4. `vortexflow/geometry.py` — builds PyLBM geometry objects.
5. `vortexflow/lbm.py` — thin adapter that configures and calls PyLBM.
6. `vortexflow/postprocess.py` — calculates engineering results.
7. `vortexflow/numerics.py` — implements course numerical algorithms.
8. `vortexflow/optimization.py` — constructs and optimizes the objective.

You do not need to understand PyLBM’s internal source code to explain this project.

## 8. Install the project

Install Python packages:

```bash
python -m pip install -r requirements.txt
```

PyLBM imports MPI even for a one-process simulation. On Ubuntu/WSL, install its runtime once:

```bash
sudo apt-get update
sudo apt-get install libopenmpi40
```

The exact package name may differ on older Ubuntu versions; `libopenmpi3` is common there.

## 9. Check the project

```bash
python verify_project.py
```

The tests cover a real short PyLBM run and all independent numerical-method functions.

## 10. Run one simulation

```bash
python run_project.py simulate configs/smoke.json
```

Results are written under:

```text
data/raw/pylbm_smoke/
```

Important files are:

| File | Meaning |
|---|---|
| `config_resolved.json` | Exact input and derived settings |
| `geometry_summary.json` | Important shape dimensions |
| `run_status.json` | Library, scheme and quality status |
| `time_series.csv` | Vortex signal and pressure history |
| `metrics.json` | Final calculated quantities |
| `field_*.npz` | Saved density and velocity fields |
| `figures/signals.png` | Vortex signal and pressure loss |
| `figures/spectrum.png` | FFT frequency spectrum |
| `figures/flow_fields.png` | Velocity and vorticity fields |

## 11. What do we receive from PyLBM?

At each requested output time, PyLBM gives three main arrays:

```text
rho = density at every fluid cell
qx  = x-direction momentum at every fluid cell
qy  = y-direction momentum at every fluid cell
```

We calculate velocity with:

```text
ux = qx / rho
uy = qy / rho
```

We calculate lattice pressure with:

```text
p = (rho - 1) / 3
```

Everything after this point is ordinary Python/NumPy analysis.

## 12. The vortex sensor signal

The earlier custom solver calculated momentum-exchange lift. That complexity has been removed.

The simplified project uses two pressure sensors in the wake:

```text
                 upper sensor •
flow ---> triangle >
                 lower sensor •
```

The vortex signal is:

```text
vortex_signal = upper_pressure - lower_pressure
```

Alternating vortices cause this value to move above and below zero. This is easy to explain and is consistent with how a vortex flowmeter can detect shedding.

## 13. Results calculated from the signal

Peak-to-peak frequency:

```text
frequency = 1 / average_time_between_peaks
```

FFT frequency:

```text
time signal -> FFT -> strongest resolved frequency
```

Signal amplitude:

```text
amplitude = (maximum - minimum) / 2
```

Signal symmetry:

```text
symmetry = abs(maximum + minimum) / abs(maximum - minimum)
```

## 14. Pressure loss and vorticity

Pressure loss is measured using an upstream and a downstream sensor:

```text
pressure_drop = upstream_pressure - downstream_pressure
```

Its time average is calculated independently with trapezoidal and Simpson integration.

Vorticity is calculated from the returned velocity field:

```text
vorticity = d(uy)/dx - d(ux)/dy
```

Alternating positive and negative vorticity behind the object indicates a vortex street.

## 15. Quality flags

Before using a case for optimization, check:

| Flag | Required value |
|---|---|
| `stable` | `true` |
| `frequency_reliable` | `true` |
| `sampling_stationary` | `true` |
| `low_mach_valid` | `true` |
| `mass_conservation_ok` | `true` |

The optimizer refuses cases that contain a failed flag.

## 16. Run the geometry optimization

```bash
python run_project.py sweep configs/course.json --radii 40 35 30 25 20
```

For each radius, the program asks PyLBM to construct and solve a new geometry. It then forms a score that:

- rewards frequency;
- rewards vortex-signal amplitude;
- penalizes pressure loss;
- penalizes signal asymmetry.

The response is studied with Lagrange interpolation, Newton interpolation, a natural cubic spline and a QR least-squares polynomial.

Golden-Section Search and Newton’s method search the QR curve for its maximum.

## 17. Verify the optimum

```bash
python run_project.py verify-optimum configs/course.json results/optimization/optimization_summary.json
```

This is essential. The curve optimum is only a prediction; a new PyLBM simulation must confirm it.

## 18. What to say during the presentation

> We use PyLBM as the tested flow-solver engine. It performs D2Q9 central-moment MRT collision, streaming and boundary treatment. Our Python implementation defines the vortex-flowmeter geometry, extracts density and momentum, converts them into pressure and velocity, measures vortex frequency and pressure loss, implements the course numerical methods, optimizes the geometry and verifies the predicted optimum.

## 19. What not to claim

Do not say:

- that this is OpenLB;
- that we wrote LBM collision and streaming ourselves;
- that the model is D3Q19 or three-dimensional;
- that its absolute pressure values reproduce the paper;
- that a fitted optimum is verified before the direct run.

Say clearly that it is a Python/PyLBM D2Q9 reduced model following the paper’s overall numerical workflow.
