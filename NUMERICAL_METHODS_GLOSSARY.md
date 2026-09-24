# Numerical Methods Glossary and Complete Workflow

This document explains the symbols, equations, numerical methods, and optimization steps used in the vortex-flowmeter project.

## 1. The main idea

The project studies how the curved-generator radius affects vortex shedding and pressure loss.

The optimization variable is:

\[
R = \text{side radius of the generator (mm)}
\]

For every candidate radius, a complete PyLBM simulation is performed. The simulation results are then used to construct an objective function and select the best radius.

```text
candidate radii
    ↓
one LBM simulation per radius
    ↓
frequency, pressure loss, signal amplitude, symmetry
    ↓
objective score J(R)
    ↓
interpolation and curve fitting
    ↓
Golden-Section/Newton search
    ↓
fresh LBM verification at predicted radius
```

## 2. Symbols and their meanings

| Symbol | Meaning | Project meaning |
|---|---|---|
| \(R\) | Design variable | Side radius being optimized |
| \(R_i\) | Candidate radius | One radius in the sweep, such as 20 mm |
| \(R^*\) | Optimum radius | Radius predicted or verified as best |
| \(x,y\) | Spatial coordinates | Lattice-grid location |
| \(t\) | Time | Physical or simulation time |
| \(\Delta t\) | Time step | Time between samples |
| \(\rho\) | Density | Fluid density calculated by LBM |
| \(\mathbf{u}\) | Velocity vector | Local fluid velocity |
| \(U\) | Reference velocity | Inlet/physical flow velocity |
| \(p\) | Pressure | Pressure reconstructed from density |
| \(D\) | Characteristic width | Reference body width used for Strouhal number |
| \(f\) | Frequency | Vortex-shedding frequency |
| \(St\) | Strouhal number | Dimensionless shedding frequency |
| \(A\) | Signal amplitude | Amplitude of the wake pressure signal |
| \(S\) | Signal asymmetry | Difference between signal halves or cycles |
| \(\Delta p\) | Pressure drop | Upstream pressure minus downstream pressure |
| \(J(R)\) | Objective function | Score used to compare radii |
| \(f_i\) | LBM population | Particle-distribution value in direction \(i\) |
| \(\mathbf{c}_i\) | Lattice velocity | Discrete D2Q9 direction |
| \(w_i\) | Lattice weight | Weight associated with direction \(i\) |
| \(c_s\) | Lattice sound speed | Sound speed used by the LBM model |
| \(\mathbf{M}\) | Moment transformation matrix | Converts populations to moments |
| \(\mathbf{S}\) | Relaxation matrix | Controls moment relaxation in MRT |

## 3. Input parameters

The main inputs are in `configs/course.json`.

| Parameter | Meaning |
|---|---|
| `nx`, `ny` | Number of lattice cells in the x and y directions |
| `steps` | Number of LBM time steps |
| `sample_start` | Time step at which transient data collection begins |
| `sample_interval` | Number of steps between saved samples |
| `reynolds_number` | Flow regime indicator, \(Re=UD/\nu\) |
| `blockage_ratio` | Generator height relative to pipe height |
| `generator_length_ratio` | Generator length relative to its height |
| `side_radius_mm` | Radius being tested or optimized |
| `fillet_radius_mm` | Radius of the downstream rounded cap |
| `incoming_angle_deg` | Flow direction angle |
| `physical_velocity_mps` | Physical inlet velocity |
| `pipe_diameter_m` | Physical pipe diameter |
| `physical_density_kgpm3` | Physical fluid density |

The nominal volumetric flow rate is calculated from the configured inlet speed and circular pipe area:

\[
Q=UA=U\frac{\pi D_{pipe}^2}{4}
\]

The corresponding nominal mass flow rate is:

\[
\dot m=\rho Q
\]

These are reported as reference flow rates. The current reduced 2-D LBM domain does not directly integrate a three-dimensional inlet flux.

## 4. Lattice Boltzmann equations

PyLBM performs the collision and streaming operations internally. The basic LBM update is:

\[
f_i(\mathbf{x}+\mathbf{c}_i\Delta t,t+\Delta t)
=
f_i(\mathbf{x},t)+\Omega_i
\]

For the central-moment MRT/Geier model, the collision operator is written conceptually as:

\[
\Omega
=
-\mathbf{M}^{-1}\mathbf{S}
\left(\mathbf{M}f-\mathbf{M}f^{eq}\right)
\]

The equilibrium distribution is:

\[
f_i^{eq}=w_i\rho\left[
1+\frac{\mathbf{c}_i\cdot\mathbf{u}}{c_s^2}
+\frac{(\mathbf{c}_i\cdot\mathbf{u})^2}{2c_s^4}
-\frac{\mathbf{u}\cdot\mathbf{u}}{2c_s^2}
\right]
\]

Macroscopic variables are recovered by:

\[
\rho=\sum_i f_i
\]

\[
\rho\mathbf{u}=\sum_i f_i\mathbf{c}_i
\]

The pressure relation is approximately:

\[
p=c_s^2\rho
\]

## 5. Time-series measurements

Two wake sensors measure pressure above and below the body. The vortex signal is:

\[
s(t;R)=p_{upper}(t;R)-p_{lower}(t;R)
\]

The pressure-loss signal is:

\[
\Delta p(t;R)=p_{upstream}(t;R)-p_{downstream}(t;R)
\]

The program stores these values in `time_series.csv`.

## 6. Frequency and FFT

Vortex shedding produces an oscillating signal. The FFT converts the signal from time domain to frequency domain:

```text
pressure versus time  →  amplitude versus frequency
```

The largest non-zero FFT peak gives the dominant frequency:

\[
f(R)=\text{dominant vortex-shedding frequency}
\]

The Strouhal number is:

\[
St(R)=\frac{f(R)D}{U}
\]

FFT is used because it identifies the dominant periodic frequency objectively, even when the time signal contains small fluctuations.

## 7. Numerical integration

The average pressure drop is calculated using the trapezoidal rule:

\[
I_{trap}\approx\sum_{j=0}^{n-2}
\frac{\Delta p_j+\Delta p_{j+1}}{2}(t_{j+1}-t_j)
\]

Simpson's rule provides an independent estimate:

\[
I_{Simp}=\frac{h}{3}
\left[y_0+y_n+4\sum y_{odd}+2\sum y_{even}\right]
\]

The integration difference is:

\[
E_I=|I_{trap}-I_{Simp}|
\]

This is an estimate of numerical integration error.

## 8. Objective function

For every radius, the project calculates:

\[
J(R)=
0.35\frac{f(R)}{f_0}
+0.20\frac{A(R)}{A_0}
-0.30\frac{\overline{\Delta p}(R)}{\overline{\Delta p}_0}
-0.15\frac{S(R)}{S_0}
\]

The terms mean:

| Term | Effect |
|---|---|
| Frequency term | Rewards a strong, measurable shedding frequency |
| Signal amplitude term | Rewards a clear wake signal |
| Pressure-loss term | Penalizes excessive pressure loss |
| Asymmetry term | Penalizes an unbalanced or irregular signal |

The baseline quantities with subscript 0 come from the first row of the sweep. The optimum is:

\[
R^*=\arg\max_R J(R)
\]

## 9. Interpolation and curve fitting

The sweep creates discrete points:

\[
(R_1,J_1),(R_2,J_2),\ldots,(R_n,J_n)
\]

The following methods estimate values between the simulated radii:

- **Lagrange interpolation:** one polynomial passes through all data points.
- **Newton interpolation:** equivalent polynomial form using divided differences.
- **Natural cubic spline:** joins local cubic polynomials smoothly.
- **QR least-squares fitting:** finds a polynomial that approximately fits the data and reduces sensitivity to noise.

These methods do not run LBM. They approximate the response curve \(\widehat{J}(R)\).

## 10. Optimization methods

Golden-Section Search maximizes \(\widehat{J}(R)\) within the tested radius interval.

Newton's optimization method searches for a stationary point:

\[
\frac{d\widehat{J}}{dR}=0
\]

Both methods produce a predicted radius. The prediction must be checked by a new LBM simulation.

## 11. Verification

If the optimizer predicts:

\[
R^*=48\text{ mm}
\]

the program runs a new simulation at exactly 48 mm:

```text
LBM(48 mm) → direct frequency, pressure, signal and objective
```

The direct result is compared with the interpolated prediction:

\[
E_{verify}=|J_{LBM}(48)-\widehat{J}(48)|
\]

Verification does not replace the sweep. It checks whether the predicted optimum behaves correctly in a fresh simulation.

## 12. Numerical errors and validity checks

The project considers:

- **Round-off error:** finite-precision computer arithmetic.
- **Truncation error:** replacing continuous equations with discrete formulas.
- **Spatial discretization error:** finite lattice resolution.
- **Temporal discretization error:** finite LBM time steps.
- **FFT resolution error:** limited by the time-window length.
- **Integration error:** difference between trapezoidal and Simpson results.
- **Sampling error:** insufficiently settled wake signal.

The code checks:

- stable simulation;
- low lattice Mach number;
- mass conservation;
- stationary sampling window;
- reliable dominant frequency.

Unreliable cases are rejected by the optimizer instead of being used to select a radius.

Because the model has open inlet and outlet boundaries, total mass in the rectangular computational domain can drift slightly as fluid enters and leaves. The measured drift is still reported, but the current acceptance limit is 2%:

\[
\frac{|M_{final}-M_{initial}|}{M_{initial}}\le0.02
\]

This is different from claiming that every cell has constant mass; it is a global domain diagnostic for an open-flow simulation.

## 13. Commands

Single fixed-geometry validation:

```bash
python run_project.py simulate configs/smoke.json
```

Actual radius sweep and optimization:

```bash
python run_project.py sweep configs/course.json --radii 20 25 30 35 40
```

Fresh verification at the predicted optimum:

```bash
python run_project.py verify-optimum \
  configs/course.json \
  results/optimization/optimization_summary.json
```

## 14. What is not used in the current project

The project does not currently need bisection, false position, Bairstow's method, Gauss elimination, Gauss-Jordan elimination, LU decomposition, power-method eigenvalues, Monte Carlo sampling, Euler, or Runge-Kutta. These methods solve different types of problems and would be artificial additions unless the project is expanded.

## 15. Exact single-case numerical sequence

For the course configuration (`steps=36000`, `sample_start=18000`, `sample_interval=5`), the solver performs the following:

1. The LBM loop advances from step 0 through step 36,000, giving 36,001 solver states.
2. Only steps 18,000 through 36,000 are sampled for the steady wake.
3. The sampled step sequence is 18,000, 18,005, ..., 36,000.
4. Therefore the number of recorded time samples is:

\[
N=\frac{36000-18000}{5}+1=3601
\]

5. The time-series file contains 3,600 equal time intervals between those 3,601 samples.
6. Trapezoidal integration uses all 3,601 points and all 3,600 intervals:

\[
I_{trap}=\sum_{j=0}^{3599}\frac{y_j+y_{j+1}}{2}\Delta t
\]

7. Simpson's rule also uses all 3,601 points, grouped into 1,800 pairs of intervals. It is valid because 3,600 is even:

\[
I_{Simp}=\frac{\Delta t}{3}
\left[y_0+y_{3600}+4\sum_{j\ odd}y_j+2\sum_{j\ even,\;j\ne0,3600}y_j\right]
\]

8. The pressure integral is divided by the observation duration (T=t_{3600}-t_0) to obtain the mean pressure drop.

## 16. Exact frequency calculations

The 3,601 equally spaced vortex-signal samples are first linearly detrended using a two-column least-squares model:

\[
s(t)\approx a_0+a_1t
\]

The fitted trend is subtracted before frequency analysis.

### FFT path

1. Apply a Hann window to the detrended 3,601-point signal.
2. Compute the real FFT (`rfft`), retaining the non-negative frequencies.
3. The frequency-bin spacing is:

\[
\Delta f=\frac{1}{T}
\]

4. Frequencies below (2/T) are ignored because they cannot contain two complete cycles in the observation window.
5. The largest remaining FFT amplitude is the FFT frequency.

### Peak-to-peak path

1. Estimate a minimum peak separation from the FFT frequency.
2. A candidate peak must be greater than its two neighboring samples.
3. Its height must exceed 15% of the detrended signal RMS.
4. Peaks closer than the estimated minimum distance are merged, retaining the larger peak.
5. At least three peaks are required.
6. Consecutive peak periods are computed; periods farther than 25% from their median are rejected.
7. At least two accepted periods are required.
8. The peak frequency is (1/\text{mean period}), and its uncertainty is the sample standard deviation of the accepted individual frequencies.

The peak frequency is marked reliable only when it agrees with the FFT frequency within the larger of (2\Delta f) or 20% of the FFT frequency.

## 17. Exact stationarity checks

The sampled series is split into two equal halves. For pressure:

\[
E_p=\frac{|\overline{p}_{second}-\overline{p}_{first}|}
{\max(|\overline{\Delta p}_{trap}|,10^{-14})}
\]

For the detrended vortex signal, the RMS of each half is calculated:

\[
E_s=\frac{|RMS_{second}-RMS_{first}|}
{\max(RMS_{first},RMS_{second},10^{-14})}
\]

The sampling window is accepted only when:

\[
E_p\le0.15\quad\text{and}\quad E_s\le0.25
\]

This is why a radius can finish its LBM run but still be rejected by the optimizer.

## 18. Exact interpolation details for a five-radius sweep

For five valid radii, the input to interpolation is five points:

\[
(R_1,J_1),\ldots,(R_5,J_5)
\]

- Lagrange uses one degree-4 polynomial through all five points.
- Newton interpolation uses the same degree-4 polynomial represented with divided differences.
- The natural cubic spline uses four cubic intervals, one between each neighboring radius pair, with zero second derivative at both endpoints.
- QR least-squares uses a degree-3 polynomial because the implementation sets:

\[
degree=\min(3,N-1)=3
\]

  It builds a five-by-four Vandermonde design matrix with columns (1,R,R^2,R^3), computes a modified Gram-Schmidt QR factorization, and solves the resulting triangular system.
- The response curves are evaluated at 401 equally spaced radii between the smallest and largest tested radius for plotting.

The 401 plotted values are not new LBM results; they are evaluations of the interpolation or fitted formulas.

## 19. Exact optimization details

Golden-Section Search operates on the QR polynomial over the interval:

\[
[R_{min},R_{max}]
\]

It repeatedly evaluates two interior points and removes the subinterval with the lower objective value. It stops when the interval width is at most (10^{-7}) mm or after 200 iterations.

Newton optimization differentiates the fitted polynomial and searches for a stationary point of:

\[
\widehat{J}'(R)=0
\]

The result is constrained to the tested radius interval. Both results are predictions from the surrogate curve.

## 20. What verification actually calculates

If the curve predicts (R^*=48) mm, verification creates a new configuration with `side_radius_mm=48.0` and repeats the full LBM process. It obtains a new direct value (J_{LBM}(48)), not an interpolated value. The prediction error is:

\[
E_{verify}=|J_{LBM}(48)-\widehat{J}(48)|
\]

The verified case must also pass the frequency, stability, low-Mach, mass-conservation, and stationarity checks.
