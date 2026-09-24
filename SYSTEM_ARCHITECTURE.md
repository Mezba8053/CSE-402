# System Architecture — PyLBM Version

## Design principle

PyLBM owns the complete LBM numerical engine. Project code treats it as a dependency and does not duplicate its algorithms.

```text
                   User command
                        |
                        v
              config.py + JSON file
                        |
          +-------------+-------------+
          |                           |
          v                           v
 geometry.py                     lbm.py adapter
 PyLBM elements              Geier scheme settings
          |                           |
          +-------------+-------------+
                        v
              PyLBM Simulation object
         collision / streaming / boundaries
                        |
                        v
                rho, qx and qy arrays
                        |
                        v
              Project sensor sampling
                        |
        +---------------+----------------+
        |                                |
        v                                v
 time_series.csv                    field_*.npz
        |                                |
        +---------------+----------------+
                        v
                 postprocess.py
                        |
     frequency / pressure / vorticity / quality
                        |
                        v
                   metrics.json
                        |
                        v
            numerics.py + optimization.py
                        |
                        v
             fitted and verified optimum
```

## Component responsibilities

| Component | Responsibility |
|---|---|
| PyLBM | Lattice, moments, equilibrium, collision, streaming, boundary algorithms and time stepping |
| `config.py` | Validate user inputs and calculate unit scales |
| `geometry.py` | Build PyLBM triangle and circle elements |
| `lbm.py` | Configure PyLBM, run it, sample its returned fields and save data |
| `pylbm_compat.py` | Restore a missing PyLBM boundary argument on Python 3.14 |
| `postprocess.py` | Calculate engineering measurements and validity flags |
| `numerics.py` | Implement course algorithms independently |
| `optimization.py` | Construct and optimize the geometry response |
| `pipeline.py` | Coordinate individual cases, sweeps and verification |
| `plots.py` | Convert saved arrays and histories into figures |

## PyLBM input contract

The adapter gives PyLBM:

- a rectangular two-dimensional domain;
- triangle and circle geometry elements;
- grid spacing;
- inlet velocity and Reynolds-number-based viscosity;
- a D2Q9 Geier central-moment scheme definition;
- Bouzidi bounce-back boundaries;
- a Neumann outlet;
- initial density and momentum.

The scheme definition is configuration required by PyLBM. The collision, streaming and boundary operations are executed inside the library.

## PyLBM output contract

The adapter reads only the conserved moments:

```text
rho  density
qx   horizontal momentum
qy   vertical momentum
```

It derives:

```text
ux = qx / rho
uy = qy / rho
p  = (rho - 1) / 3
```

The adapter then converts lattice velocity and pressure to configured physical units.

## Geometry design

The body is the combination of four PyLBM elements:

1. one solid triangle;
2. one fluid circle cutting the upper concave side;
3. one fluid circle cutting the lower concave side;
4. one solid circle making the rounded downstream cap.

Changing `side_radius_mm` changes the two cutting circles. Changing `fillet_radius_mm` changes the downstream cap.

## Sensor design

Three pressure measurements are derived from the density field:

- upstream pressure;
- downstream pressure;
- upper-minus-lower wake pressure.

The last value is the periodic vortex signal. This replaces the earlier manual momentum-exchange force implementation.

## Numerical-method layer

The numerical-method layer remains independent of PyLBM:

```text
time integration:       trapezoidal and Simpson
frequency:              peak intervals and FFT
spatial derivative:     central-difference vorticity
interpolation:           Lagrange, Newton, natural spline
response fitting:       modified Gram-Schmidt QR
optimization:           Golden-Section and Newton
```

This independence makes the course algorithms visible and unit-testable even though a library supplies CFD.

## Quality gates

Optimization accepts a simulation only when its available flags pass:

- solver remained finite;
- peak and FFT frequencies agree;
- both halves of the sample window are stationary;
- maximum lattice Mach is at most 0.30;
- relative mass change is at most 1%.

## Model boundary

The PyLBM version is a 2D far-field bluff-body model. The paper is a 3D OpenLB methane-pipe model. Comparisons should emphasize trends rather than identical absolute values.
