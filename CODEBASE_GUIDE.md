# Codebase Guide — Simplified PyLBM Edition

## The only idea you need first

`vortexflow/lbm.py` does not solve LBM equations itself. It prepares a dictionary, creates `pylbm.Simulation`, calls `one_time_step()`, and reads density and momentum.

## Execution order

When you run:

```bash
python run_project.py simulate configs/smoke.json
```

the calls are:

```text
run_project.py
  -> vortexflow.cli.main
  -> config.load_config
  -> pipeline.run_case
  -> lbm.run_simulation
       -> geometry.build_pylbm_elements
       -> pylbm.Simulation
       -> solution.one_time_step repeatedly
       -> save time series and fields
  -> postprocess.analyze_case
  -> plots.plot_case
```

## `config.py`

`SimulationConfig` holds all settings. Derived properties calculate grid spacing, physical time step, viscosity and relaxation information. Validation catches unsafe or meaningless inputs before PyLBM starts.

## `geometry.py`

This file uses public PyLBM geometry classes. It does not decide which cells collide or stream; PyLBM constructs that domain.

`build_pylbm_elements()` returns a list containing a triangle and three circles. `_arc_center()` is ordinary circle geometry used to position the concave cuts.

## `lbm.py`

Important functions are:

- `_geier_scheme()` — builds the PyLBM input dictionary;
- `_fields()` — reads `rho`, `qx` and `qy` from PyLBM;
- `_probe_values()` — samples pressure sensors;
- `run_simulation()` — controls output timing and files.

You are not expected to explain every symbolic expression in `_geier_scheme()`. It is the scheme configuration adapted from PyLBM’s official `Karman_vortex_street` example. The actual solver routines remain inside PyLBM.

## `pylbm_compat.py`

Python 3.14 changed how dynamic local variables are exposed. PyLBM 0.11 expects the older behavior when preparing boundary arrays. The compatibility file adds those array references explicitly.

It does not alter equilibrium, collision, streaming or any physical result.

## `postprocess.py`

This reads saved CSV/NPZ files and calculates:

- frequency from peaks;
- frequency from FFT;
- pressure averages;
- signal amplitude and symmetry;
- vorticity;
- stationarity and agreement flags.

Analysis can therefore be repeated without rerunning PyLBM.

## `numerics.py`

These are our course implementations:

- `trapezoidal`;
- `simpson_uniform`;
- `lagrange_evaluate`;
- `newton_coefficients` and `newton_evaluate`;
- `natural_cubic_spline`;
- `qr_least_squares`;
- `golden_section_maximize`;
- `newton_maximize_polynomial`.

## `optimization.py`

The objective rewards frequency and signal amplitude and penalizes pressure loss and asymmetry. It fits the result against side radius, runs two optimizers and saves their histories.

## `pipeline.py`

- `run_case()` runs and analyzes one case.
- `run_sweep()` repeats cases at several side radii and optimizes them.
- `verify_optimum()` performs the required new simulation at the fitted optimum.

## Safe changes

- Change runtime or resolution in JSON configuration files.
- Change geometry only in `geometry.py`.
- Add measurements in `postprocess.py` and `pipeline.METRIC_COLUMNS`.
- Change objective weights in `optimization.DEFAULT_WEIGHTS`.
- Do not edit installed PyLBM files under `.python_packages`.

## Debugging

If PyLBM cannot import `mpi4py`, install the system OpenMPI runtime.

If a run finishes but is unusable, inspect `run_status.json`, `metrics.json` and the three figures. Never replace a failed measurement with a paper value.

If a geometry sweep stops, fix the failed case rather than allowing the optimizer to use invalid data.
