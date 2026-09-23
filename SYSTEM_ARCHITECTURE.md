# System Architecture

## 1. Purpose and design boundary

The system is a fully Python numerical simulation and optimization pipeline for a triangular vortex flowmeter. It converts a version-controlled JSON case into flow histories, fields, engineering metrics, a fitted geometry response and a proposed optimum.

The production boundary is explicit:

- implemented: 2D D2Q9 MRT-LBM with a Smagorinsky correction;
- reference only: the paper’s 3D D3Q19 OpenLB implementation;
- not claimed: exact high-Reynolds-number quantitative reproduction or experimental validation.

## 2. Logical architecture

```text
                         User / CLI
                             |
                             v
                  +---------------------+
                  | Configuration layer |
                  | config.py + JSON    |
                  +----------+----------+
                             |
                             v
             +---------------+----------------+
             |                                |
             v                                v
    +------------------+             +------------------+
    | Geometry service |             | Unit/scaling     |
    | geometry.py      |             | config properties|
    +---------+--------+             +---------+--------+
              |                                |
              +---------------+----------------+
                              v
                  +-----------------------+
                  | Transient flow solver |
                  | lbm.py                |
                  +-----------+-----------+
                              |
                  +-----------+-----------+
                  |                       |
                  v                       v
          time_series.csv             field_*.npz
                  |                       |
                  +-----------+-----------+
                              v
                  +-----------------------+
                  | Analysis layer        |
                  | postprocess.py        |
                  +-----------+-----------+
                              |
               metrics.json, spectrum, vorticity
                              |
                  +-----------+-----------+
                  |                       |
                  v                       v
            plots.py              repeated cases
                                          |
                                          v
                               geometry_metrics.csv
                                          |
                                          v
                             +-----------------------+
                             | Numerical optimization|
                             | numerics.py +         |
                             | optimization.py       |
                             +-----------+-----------+
                                         |
                                         v
                              optimum + histories
                                         |
                                         v
                              direct verification run
```

## 3. Components and responsibilities

| Component | Responsibility | Must not do |
|---|---|---|
| `config.py` | Load/validate JSON and derive lattice/physical scales | Run simulations |
| `geometry.py` | Build Boolean wall and generator masks | Modify solver state |
| `lbm.py` | Equilibrium, MRT collision, LES correction, streaming, boundaries and sampling | Fit response curves |
| `postprocess.py` | Signal validation, frequency, integration and vorticity | Hide unreliable results |
| `numerics.py` | Standalone course algorithms | Read project files |
| `optimization.py` | Objective, response models and bounded searches | Claim surrogate optimum is verified |
| `plots.py` | Deterministic PNG generation with a headless backend | Change numerical results |
| `pipeline.py` | Orchestrate cases and sweeps | Contain solver mathematics |
| `cli.py` | Parse commands and display summaries | Duplicate business logic |

This separation lets numerical algorithms be unit-tested without running computational fluid dynamics, while the short solver test checks integration of the core flow components.

## 4. Solver data model

The primary state is `f[q, y, x]`, where `q = 9`. Derived two-dimensional arrays are:

- `rho[y, x]`: lattice density;
- `ux[y, x]`, `uy[y, x]`: lattice velocity;
- `solid[y, x]`: walls plus generator;
- `body[y, x]`: generator only, used for force integration.

The direction ordering is:

```text
0 rest
1 east, 2 north, 3 west, 4 south
5 northeast, 6 northwest, 7 southwest, 8 southeast
```

The moment transform `M` and its inverse are constants. MRT relaxation is local only for the two shear moments; other non-conserved moments use fixed relaxation rates.

## 5. Per-step execution sequence

```text
1. Recover rho and velocity from f
2. Check finite values and lattice-speed safety limit
3. Build equilibrium distributions
4. Calculate non-equilibrium stress
5. Calculate local Smagorinsky effective tau
6. Transform f and equilibrium to moment space
7. Relax non-conserved moments
8. Transform back to distribution space
9. Calculate generator force by momentum exchange
10. Stream distributions and bounce back solid links
11. Apply ramped Zou-He inlet and density outlet
12. Sample forces/pressure if the configured interval matches
13. Save a field snapshot if requested
```

The transverse perturbation is active only while the inlet ramps. It breaks exact mathematical symmetry but does not impose a frequency on the mature wake.

## 6. Configuration contract

Every simulation is defined by a JSON object accepted by `SimulationConfig`. Important groups are:

| Group | Fields |
|---|---|
| Identity/output | `case_id`, `output_root`, `save_snapshots` |
| Grid/time | `nx`, `ny`, `steps`, `sample_start`, `sample_interval`, `snapshot_interval` |
| LBM | `inlet_lattice_velocity`, `reynolds_number`, `smagorinsky_constant`, `inlet_perturbation_fraction` |
| Geometry | `blockage_ratio`, `generator_x_ratio`, `generator_length_ratio`, `incoming_angle_deg`, `side_radius_mm`, `fillet_radius_mm` |
| Physical scale | `pipe_diameter_m`, `physical_velocity_mps`, `physical_density_kgpm3`, `characteristic_width_m` |

The resolved configuration also records `dx`, `dt`, generator dimensions, relaxation time and model name. A case must be reproducible using this resolved file.

## 7. File contracts

### Time-series contract

`time_series.csv` contains one row per retained sample:

```text
step,time_s,lift_npm,drag_npm,p_upstream_pa,p_downstream_pa,
pressure_drop_pa,max_lattice_velocity
```

`lift_npm` and `drag_npm` are forces per unit depth because the solver is two-dimensional.

### Field contract

Each compressed snapshot contains:

```text
rho, ux, uy, solid, body, step, time_s
```

Velocities in snapshot files are physical `m/s`; density is the lattice density used for the isothermal pressure relation.

### Metrics contract

Case metrics include:

- geometry and Reynolds number;
- peak and FFT frequencies, uncertainty, resolution and reliability;
- Strouhal number, lift amplitude and symmetry deviation;
- trapezoidal and Simpson pressure means;
- maximum absolute vorticity;
- stability, mass-change and maximum-speed diagnostics.

Optimization must reject missing/non-finite required values and unreliable simulation inputs.

## 8. Analysis design

### Frequency

Frequency extraction has two independent paths:

```text
lift -> linear detrend -> Hann FFT -> dominant resolvable frequency
                    |
                    `-> frequency-informed peak spacing -> mean period
```

The observation window must include at least two full cycles for an FFT candidate and at least three peaks for a peak estimate. Peak-period scatter must be at most 20%. The final reliability flag also requires agreement with FFT within 20% or two FFT bins.

### Pressure integration

The same sampled pressure-drop history is integrated by trapezoidal and Simpson rules. Simpson’s result is only defined for uniform sampling with an even number of intervals. Configuration files should choose start, stop and interval values accordingly.

### Spatial derivative

Vorticity is evaluated from physical velocity snapshots using second-order central differences in the interior through `numpy.gradient`. Solid cells are masked before statistics and plotting.

## 9. Optimization architecture

Input cases share identical flow/grid settings and differ only in the selected geometry parameter. The baseline is the first row of the sweep table. Each physical response is normalized before the weighted objective is formed.

Four response representations are saved:

- Lagrange interpolating polynomial;
- Newton divided-difference interpolating polynomial;
- natural cubic spline;
- degree-at-most-three least-squares polynomial solved by modified Gram-Schmidt QR.

Golden-Section Search and Newton’s method maximize the QR polynomial on the simulated interval. Both histories are retained. The chosen radius is not considered final until a fresh flow simulation validates it.

## 10. Failure behavior

The design fails explicitly when:

- a configuration is outside supported bounds;
- relaxation time is too close to `0.5`;
- distributions or velocity become non-finite;
- maximum lattice velocity exceeds `0.5`;
- the time series is too short;
- the signal is at round-off scale;
- peak periods are inconsistent;
- an optimization input is missing/non-finite;
- fewer than four sweep cases are supplied.

A stable solver status alone does not establish a valid shedding frequency.

## 11. Validation gates

| Gate | Acceptance evidence |
|---|---|
| Unit algorithms | `verify_project.py` passes |
| Solver smoke | finite state, speed below limit, result files produced |
| Wake physics | alternating vorticity downstream and non-roundoff lift |
| Signal quality | multiple mature cycles and peak/FFT agreement |
| Integration | trapezoidal/Simpson difference is reported and small enough for the study |
| Grid sensitivity | trends persist on at least three grids |
| Sweep consistency | identical non-geometry settings and reliable metrics |
| Optimizer agreement | Golden and Newton predictions are compared |
| Final verification | predicted optimum is simulated directly |

## 12. Extension rules

- Add a new configuration field in `SimulationConfig`, validation, both example JSON files and this document.
- Add a metric in `analyze_case`, `METRIC_COLUMNS`, plotting/reporting and tests.
- Add a geometry parameter only through `geometry.py`; do not place shape logic in the time loop.
- Keep numerical methods pure and file-independent in `numerics.py`.
- Never replace a failed measurement with a plausible constant or a paper value.
- Do not compare 2D reduced-model absolute values with 3D paper values without stating the modelling difference.

For a file-by-file reading order, see [CODEBASE_GUIDE.md](CODEBASE_GUIDE.md).
