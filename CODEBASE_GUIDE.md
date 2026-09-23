# Codebase Guide

This guide explains where to start, how control and data move through the program, and where to make common changes.

## 1. Recommended reading order

1. `configs/smoke.json` — one complete input case.
2. `run_project.py` — executable entry point and local dependency discovery.
3. `vortexflow/cli.py` — available commands.
4. `vortexflow/pipeline.py` — high-level case and sweep workflows.
5. `vortexflow/config.py` — input contract and physical/lattice scaling.
6. `vortexflow/geometry.py` — grid masks.
7. `vortexflow/lbm.py` — flow solver.
8. `vortexflow/postprocess.py` — quantities extracted from a run.
9. `vortexflow/numerics.py` — course numerical methods.
10. `vortexflow/optimization.py` — response construction and optimization.
11. `vortexflow/plots.py` — figures.
12. `tests/` — executable examples of expected behavior.

## 2. Entrypoints

`run_project.py` inserts the optional repository-local `.python_packages` directory into `sys.path`, then calls `vortexflow.cli.main()`. A normal virtual environment works as well; the local directory is only a convenience for restricted lab machines.

`verify_project.py` discovers every `unittest` file under `tests/` and exits with a nonzero code if any check fails.

The CLI routes commands as follows:

```text
simulate -> load_config -> pipeline.run_case
analyze  -> postprocess.analyze_case
sweep    -> pipeline.run_sweep -> optimization.optimize_metrics
optimize -> optimization.optimize_metrics + optimization plot
```

## 3. One simulation from start to finish

### Configuration

`load_config()` maps JSON keys to `SimulationConfig`. Its properties derive body size, physical `dx` and `dt`, lattice viscosity, base relaxation time and output directory.

`validate()` is the first safety boundary. When a new run fails before step zero, fix the configuration rather than weakening the check without a numerical reason.

### Geometry

`generator_mask()` creates a Boolean array for the bluff body. The upstream face, body length, height, side curvature and downstream circular cap are built from configured parameters. The paper does not provide full CAD coordinates, so this shape is explicitly an engineering approximation.

`build_solid_mask()` adds the upper and lower channel walls. It returns both `solid` and `body`: the solver reflects populations from all solids but integrates force only on the generator.

### Solver initialization and collision

`run_simulation()` starts with lattice density one and zero velocity. `equilibrium()` constructs the nine initial distributions. The initial fluid mass is saved for a final diagnostic.

`macroscopic()` recovers density and velocity. `_effective_tau()` calculates a local Smagorinsky correction from non-equilibrium stress. `_collide_mrt()` then:

1. transforms distributions with `M`;
2. relaxes each moment with its assigned rate;
3. transforms back with `M_INV`.

Density and momentum moments are conserved. The two shear moments use the local effective rate.

### Force, streaming and boundaries

`_momentum_exchange()` sums momentum reflected from fluid-to-generator links. `_stream_with_bounce_back()` streams non-blocked populations and reflects blocked populations at their source fluid nodes.

`_apply_boundaries()` implements a Zou-He velocity inlet and constant-density outlet. The main loop ramps the inlet smoothly and adds a small transverse pulse only during that ramp to seed the antisymmetric wake.

### Sampling

After `sample_start`, the solver writes every `sample_interval` steps. Pressure probes are positioned upstream and downstream using the body size. Snapshots are less frequent because fields are much larger than scalar histories.

The solver stops with an error if any value becomes non-finite or lattice speed exceeds the conservative project limit.

## 4. Post-processing

`analyze_case()` reads only saved outputs, so analysis can be repeated without rerunning CFD.

Signal processing:

1. remove a least-squares linear trend;
2. apply a Hann window for FFT;
3. ignore frequencies too low to finish two cycles;
4. use the FFT candidate to prevent noise peaks from being counted too close together;
5. reject inconsistent periods;
6. compare peak and FFT estimates.

Pressure is integrated over physical time. Vorticity is computed from the most recent physical-velocity field. The resulting `metrics.json` is the authoritative summary of one case.

When `frequency_reliable` is false, inspect `signals.png`, `spectrum.png`, the sample window and the wake field. Do not manually replace the value.

## 5. Course numerical methods

`numerics.py` intentionally implements the algorithms instead of wrapping SciPy:

| Function | Algorithm |
|---|---|
| `trapezoidal` | composite trapezoidal integration |
| `simpson_uniform` | composite Simpson 1/3 integration |
| `lagrange_evaluate` | direct Lagrange basis |
| `newton_coefficients` | divided-difference table in-place |
| `newton_evaluate` | nested Newton evaluation |
| `natural_cubic_spline` | natural spline linear system |
| `qr_least_squares` | modified Gram-Schmidt QR and solve |
| `golden_section_maximize` | bounded derivative-free search |
| `newton_maximize_polynomial` | derivative iteration with bounds |

NumPy is used for arrays and linear algebra, but the course algorithms and their iteration logic remain visible.

## 6. Geometry sweep and objective

`run_sweep()` clones the base configuration with `dataclasses.replace()`, changes `case_id` and `side_radius_mm`, and runs each case independently. It collects metrics in `data/processed/geometry_metrics.csv`.

`optimize_metrics()` checks the data, forms the weighted normalized objective, produces all interpolation curves, fits the QR polynomial and runs both optimizers. Coefficients are stored in increasing-power order:

```text
J(R) = a[0] + a[1] R + a[2] R^2 + ...
```

The baseline for normalization is the first input row. Keep the intended baseline first when preparing a table manually.

## 7. Common changes

### Change resolution or duration

Edit a copied JSON file. Increase `ny` to resolve the channel and generator; increase `nx` enough to retain a downstream wake. When changing `ny`, check the derived relaxation time in `config_resolved.json`.

### Add a geometry variable

Add it to `SimulationConfig`, validate it, then use it only in `generator_mask()`. Add its value to metrics and the sweep table if it must be analyzed.

### Add an output metric

Calculate it in `analyze_case()`, add it to `pipeline.METRIC_COLUMNS`, add a focused unit test and update both documentation contracts.

### Change objective weights

Edit `DEFAULT_WEIGHTS` or pass a weight dictionary to `optimize_metrics()` from Python. Record the weights in the report; a different weighting represents a different engineering preference.

### Use a different baseline

Put that case first in the metrics CSV before optimization. A future extension could expose an explicit baseline case ID in the CLI.

## 8. Debugging checklist

If a simulation is unstable:

- reduce `inlet_lattice_velocity`;
- lower Reynolds number for the course model;
- increase generator resolution;
- verify `base_relaxation_time` is not too close to `0.5`;
- check that the outlet is sufficiently far downstream.

If no vortex frequency is accepted:

- confirm lift is larger than round-off scale;
- start sampling after inlet ramp and wake development;
- extend the run to include several cycles;
- inspect alternating positive/negative wake vorticity;
- compare FFT resolution with the expected frequency;
- do not increase the startup perturbation merely to create a desired answer.

If pressure estimates disagree:

- make the number of sampled intervals even for Simpson’s rule;
- extend the averaging window;
- inspect pressure history for remaining startup drift;
- move probes only with a documented geometric reason.

## 9. Test coverage

The suite checks equilibrium/macroscopic consistency, geometry masks, a finite short LBM run, numerical integration, all three interpolation approaches, QR polynomial recovery, both optimizers, known sinusoidal frequency, rejection of round-off pseudo-signals and vorticity of a known rotating field.

These tests establish code behavior, not physical validation of every flow regime. Grid/time convergence and comparison against literature remain required report tasks.

## 10. Scientific interpretation limits

- D2Q9 is a 2D reduced model; the paper uses D3Q19.
- The geometry is reconstructed from reported parameters, not original CAD.
- Course Reynolds numbers are lower than the methane pipeline values for numerical practicality.
- LES in two dimensions is a numerical closure demonstration, not the same turbulence physics as 3D LES.
- Absolute force and pressure values should not be presented as experimentally validated.
- Trends become defensible only after sensitivity and repeatability studies.

The safest presentation is: “We implemented the paper’s numerical workflow in a transparent Python reduced-order model, then added the interpolation, QR fitting and optimization methods required by the final proposal.”
