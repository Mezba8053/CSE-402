# Validation Report — PyLBM Refactor

Date: 23 September 2026

## What was validated

The former hand-written LBM implementation was removed. The verified flow path is now:

```text
project configuration -> PyLBM 0.11.0 -> returned rho/qx/qy fields
-> project sensors and analysis -> figures and metrics
```

## Automated tests

Command:

```bash
python verify_project.py
```

Result: `Ran 10 tests ... OK`. The suite contains a real short PyLBM simulation plus independent tests for integration, interpolation, QR fitting, optimization, frequency detection and vorticity.

## Full smoke simulation

Command:

```bash
python run_project.py simulate configs/smoke.json
```

Solver recorded in `run_status.json`:

```text
PyLBM 0.11.0
D2Q9 Geier central-moment MRT
NumPy backend
```

Configuration:

| Setting | Value |
|---|---:|
| Grid | 140 x 42 |
| Steps | 18,000 |
| Retained window | 10,000–18,000 |
| Reynolds number | 150 |
| Inlet lattice velocity | 0.035 |
| Side radius | 40 mm |
| End fillet | 4 mm |

Results:

| Metric | Value |
|---|---:|
| Peak-to-peak frequency | 101.4011 Hz |
| Peak-frequency standard deviation | 2.1988 Hz |
| FFT frequency | 110.2470 Hz |
| FFT resolution | 15.7594 Hz |
| Peak/FFT difference | 8.02% |
| Detected peaks | 7 |
| Strouhal number | 0.27808 |
| Vortex-signal amplitude | 2747.21 Pa |
| Mean pressure drop, trapezoidal | 1985.5386 Pa |
| Mean pressure drop, Simpson | 1985.5380 Pa |
| Integration difference | 0.00062 Pa |
| Maximum lattice Mach | 0.14893 |
| Relative mass change | 0.729% |

All quality gates passed:

```text
stable = true
frequency_reliable = true
sampling_stationary = true
low_mach_valid = true
mass_conservation_ok = true
```

Visual inspection confirmed an alternating positive/negative vorticity street behind the generator and a periodic upper-minus-lower pressure signal.

## Interpretation

These results validate the software workflow and the reduced PyLBM model. Absolute pressure values are not expected to equal the paper because the domain, dimensionality, Reynolds number, turbulence treatment and sensor definition differ.

A final engineering conclusion still requires the complete radius sweep, multiple grid resolutions and direct simulation of the fitted optimum.
