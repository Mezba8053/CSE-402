# Python Numerical Simulation and Optimization of a Triangular Vortex Flowmeter

This repository is an executable, course-scale numerical project based on the final proposal and the paper:

> Q. Sun et al., “Shape Optimization of the Triangular Vortex Flowmeter Based on the LBM Method,” *Symmetry*, 17(4), 534, 2025.

The project runs a transient lattice-Boltzmann flow simulation, extracts vortex-flowmeter quantities, applies the numerical methods named in the proposal, and supports a side-radius optimization study. All source code is Python.

## Important model scope

The paper uses OpenLB with a three-dimensional D3Q19 MRT-LBM and LES model. OpenLB itself is a C++ framework; using its normal solver would violate the requirement that all project code be Python. Therefore this repository implements a reduced, fully Python model:

| Item | Paper | This implementation |
|---|---|---|
| Language/framework | C++ OpenLB | Python and NumPy |
| Lattice | 3D D3Q19 | 2D D2Q9 |
| Collision | MRT | MRT |
| Turbulence treatment | LES | local Smagorinsky relaxation correction |
| Geometry | 3D pipe and generator | 2D channel cross-section and generator |
| Frequency | FFT | peak-to-peak plus FFT cross-check |
| Optimization | discrete comparison | interpolation, QR fitting, Golden Section and Newton |

This is suitable for demonstrating numerical simulation and numerical-method implementation. It is not a quantitatively exact reproduction of the paper’s high-Reynolds-number 3D OpenLB result. If the course explicitly requires OpenLB, the “all code in Python” constraint must be relaxed and a separate C++ OpenLB solver must be used.

## Current implementation status

- D2Q9 MRT collision and streaming: implemented
- Smagorinsky effective-relaxation correction: implemented
- Parameterized triangular/concave generator: implemented
- No-slip walls and obstacle bounce-back: implemented
- Zou-He velocity inlet and density outlet: implemented
- Startup-only wake perturbation: implemented
- Momentum-exchange lift and drag: implemented
- Pressure probes and physical-unit conversion: implemented
- Vorticity and Strouhal number: implemented
- Peak-to-peak and FFT frequency extraction: implemented
- Frequency reliability checks: implemented
- Trapezoidal and Simpson integration: implemented
- Lagrange, Newton and cubic-spline interpolation: implemented
- Polynomial least squares by explicit QR: implemented
- Golden-Section and Newton optimization: implemented
- Automated tests and plots: implemented

## Quick start

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
python verify_project.py
python run_project.py simulate configs/smoke.json
```

On a Windows machine where the project is run through WSL:

```bash
python3 -m pip install -r requirements.txt
MPLCONFIGDIR=/tmp/vortexflow-mpl python3 run_project.py simulate configs/smoke.json
```

The smoke case is smaller than a production study, but it still runs the complete path: configuration, geometry, LBM, CSV output, signal analysis, vorticity, metrics and figures.

## Commands

Run one case and post-process it:

```bash
python run_project.py simulate configs/smoke.json
```

Reanalyze a completed case without rerunning the solver:

```bash
python run_project.py analyze data/raw/smoke
```

Run a geometry sweep and optimize the resulting response:

```bash
python run_project.py sweep configs/course.json --radii 40 35 30 25 20
```

Optimize an existing validated table:

```bash
python run_project.py optimize data/processed/geometry_metrics.csv
```

Run all verification tests:

```bash
python verify_project.py
```

## End-to-end workflow

```text
JSON configuration
      |
      v
2D generator + channel mask
      |
      v
D2Q9 MRT-LBM + Smagorinsky correction
      |
      +--> field snapshots (.npz)
      |
      `--> lift, drag and pressure history (.csv)
                   |
                   v
       transient-window post-processing
                   |
       +-----------+------------+
       |           |            |
       v           v            v
   frequency   pressure mean  vorticity
   peak + FFT  trap + Simpson central differences
       |           |            |
       +-----------+------------+
                   |
                   v
              metrics.json
                   |
       repeated for each side radius
                   |
                   v
    interpolation + QR response fitting
                   |
                   v
       Golden-Section + Newton search
                   |
                   v
       direct optimum verification run
```

Do not accept an optimum until it has been run as a new flow case. An optimum of a fitted polynomial is only a prediction.

## What the paper does numerically

The paper treats methane flow through a DN100 pipeline containing a triangular vortex generator. At the numerical level it:

1. Builds a 3D computational pipe and bluff-body geometry.
2. Uses OpenLB’s D3Q19 lattice, MRT collision and LES treatment.
3. Applies a velocity inlet, pressure outlet and no-slip solid boundaries.
4. Performs a mesh-independence study and refines the generator/wake region.
5. Runs transient simulations until a periodic vortex street is established.
6. Compares candidate incoming angles, end fillets and concave-side radii.
7. Examines vorticity, lift history, lift symmetry, pressure loss and vortex frequency.
8. Uses an FFT of the lift signal to identify the shedding frequency.
9. Selects the best tested geometry rather than performing a continuous mathematical optimization.

The reported selected geometry is an incoming angle of 180 degrees, a concave-side radius of 25 mm and a downstream fillet radius of 4 mm. Selected paper results are:

| Velocity | Quantity | Original | Optimized |
|---:|---|---:|---:|
| 10.21 m/s | Pressure drop | 404.54 Pa | 334.93 Pa |
| 10.21 m/s | Frequency | 131.95 Hz | 135.55 Hz |
| 2.94 m/s | Pressure drop | 79.12 Pa | 36.51 Pa |
| 2.94 m/s | Frequency | 32.66 Hz | 37.18 Hz |

These values are literature references, not automatic pass/fail targets for the reduced 2D model.

## Numerical model used here

For each D2Q9 direction `i`, the solver evolves a distribution `f_i`:

```text
f_i(x + c_i, t + 1) = f_i(x,t) + collision_i
```

Macroscopic quantities are recovered using:

```text
rho   = sum_i(f_i)
rho*u = sum_i(f_i*c_i)
p     = c_s^2 (rho - rho_reference),  c_s^2 = 1/3
```

The collision is performed in moment space:

```text
m* = m - S (m - m_equilibrium)
f* = inverse(M) m*
```

Conserved density and momentum moments have zero relaxation. Shear moments use a local effective relaxation time corrected by the non-equilibrium stress and Smagorinsky constant. Populations are then streamed. Links entering a solid cell are reflected in the opposite direction.

The flow starts from rest and the inlet is smoothly ramped. A small transverse pulse exists only during startup. This is necessary because a perfectly symmetric discrete geometry and exactly symmetric initial conditions can remain on the unstable symmetric solution and never select an alternating shedding side.

### Physical/lattice conversion

The channel height represents the physical pipe diameter:

```text
dx = pipe_diameter / (ny - 2)
dt = inlet_lattice_velocity * dx / physical_velocity
```

Pressure and force per unit depth are scaled using the configured physical density. Because the model is 2D, force is reported as `N/m`, not total 3D force.

## Calculated quantities

Vorticity:

```text
omega_z = d(uy)/dx - d(ux)/dy
```

Strouhal number:

```text
St = frequency * characteristic_width / inlet_velocity
```

Lift symmetry deviation:

```text
eta = abs(Lmax + Lmin) / abs(Lmax - Lmin)
```

Time-averaged pressure loss:

```text
mean_delta_p = integral(delta_p(t) dt) / observation_duration
```

Both composite trapezoidal and composite Simpson rules are evaluated when the sampling grid satisfies Simpson’s requirements.

### Frequency validity

The post-processor removes a linear trend, finds the dominant FFT frequency and uses it to set a safe minimum distance between time-domain peaks. A peak result is accepted only if its period scatter is small and it agrees with FFT within 20% or two FFT bins. Important fields in `metrics.json` are:

- `frequency_hz`: peak-to-peak estimate;
- `frequency_std_hz`: scatter of accepted individual periods;
- `fft_frequency_hz`: independent FFT estimate;
- `fft_resolution_hz`: `1 / observation_duration`;
- `frequency_reliable`: agreement flag.

A geometry sweep should only be optimized when every input case has `frequency_reliable = true`, is stable and contains enough mature cycles.

## Optimization added by the proposal

For every tested side radius `R`, the pipeline forms a normalized score:

```text
J(R) = 0.35*(f/f0)
     + 0.20*(lift_amplitude/lift0)
     - 0.30*(pressure_drop/pressure0)
     - 0.15*(symmetry/symmetry0)
```

The weights are explicit defaults in `vortexflow/optimization.py` and may be changed. The first row of the metrics table is the baseline used for normalization.

The project then:

1. constructs Lagrange, Newton divided-difference and natural cubic-spline interpolants;
2. fits a degree-at-most-three response using modified Gram-Schmidt QR;
3. maximizes the QR polynomial using Golden-Section Search;
4. applies Newton’s method using analytical polynomial derivatives;
5. saves both iteration histories;
6. requires a direct LBM run at the proposed optimum.

## Outputs

For a case named `smoke`, the solver writes:

```text
data/raw/smoke/
|-- config_resolved.json
|-- run_status.json
|-- time_series.csv
|-- metrics.json
|-- spectrum.npz
|-- vorticity.npy
|-- field_*.npz
`-- figures/
    |-- signals.png
    |-- spectrum.png
    `-- flow_fields.png
```

Sweep and optimizer outputs are written to:

```text
data/processed/geometry_metrics.csv
results/optimization/
|-- response_curves.csv
|-- optimization_summary.json
|-- golden_history.json
|-- newton_history.json
`-- optimization.png
```

Generated result directories are ignored by Git because they can be recreated from the configurations.

## Repository map

```text
configs/                 reproducible run settings
tests/                   unit and short solver tests
vortexflow/config.py     configuration, validation and unit scales
vortexflow/geometry.py   channel and generator masks
vortexflow/lbm.py        D2Q9 MRT-LBM-LES solver
vortexflow/postprocess.py frequency, pressure and vorticity analysis
vortexflow/numerics.py   course numerical-method algorithms
vortexflow/optimization.py objective, fitting and optimizers
vortexflow/plots.py      non-interactive result figures
vortexflow/pipeline.py   case and sweep orchestration
vortexflow/cli.py        command-line interface
run_project.py           user entry point
verify_project.py        complete test entry point
CODEBASE_GUIDE.md        detailed code-reading guide
SYSTEM_ARCHITECTURE.md   component and data-contract design
```

## Validation checklist for the report

- Run `python verify_project.py` and record the result.
- Report grid dimensions, Reynolds number, `tau`, maximum lattice velocity and mass change.
- Show that the wake alternates in the vorticity plot.
- Exclude startup data and include several mature cycles.
- Report peak and FFT frequencies with FFT resolution.
- Compare trapezoidal and Simpson pressure means.
- Repeat the baseline on at least two additional grids.
- Run every sweep point with the same numerical settings.
- Simulate the predicted optimum directly.
- Compare trends with the paper, while clearly labeling the 2D/reduced-model limitation.

See [CODEBASE_GUIDE.md](CODEBASE_GUIDE.md) for the execution path and safe extension points, and [SYSTEM_ARCHITECTURE.md](SYSTEM_ARCHITECTURE.md) for the formal design.
