# TART-GNSS Stage 5 Implementation Specification
## TLE-Independent Statistical Trajectory Tracking and Bayesian State Estimation

**Repository:** `git@github.com:ifornax/tart-gnss.git`  
**Base branch:** `main`  
**Stage status before this work:** Stages 1–4 implemented  
**This specification implements:** Stage 5 only  
**Canonical Stage-5 CLI:** `tart-track`

---

# 0. Permanent CLI Naming Convention

Every public command must begin with `tart-`.

Current canonical commands:

```text
tart-stations-update
tart-stations-list
tart-common-gnss
tart-retrieve
tart-localise
tart-triangulate
```

Stage 5 must add:

```text
tart-track
```

Do not introduce non-canonical aliases such as `track-satellite`, `trajectory-fit`, or `gnss-track`. Update documentation, tests, entry points and examples consistently. Historical provenance files must not be rewritten.

---

# 1. Objective

Implement **Stage 5: statistical trajectory inference from the time-resolved 3-D positions and covariances produced by Stage 4**.

Stage 4 provides independent position measurements

\[
\hat{\mathbf r}_k
\]

at epochs

\[
t_k
\]

with measurement covariance

\[
R_k.
\]

Stage 5 must infer a continuous satellite state

\[
\mathbf x_k =
\begin{bmatrix}
\mathbf r_k \\
\mathbf v_k
\end{bmatrix}
\]

with position and velocity plus a full state covariance.

The primary outputs are:

1. filtered position and velocity estimates;
2. smoothed position and velocity estimates;
3. full state covariance;
4. innovation statistics;
5. uncertainty-calibrated quality flags;
6. diagnostic plots;
7. a clean Stage-6 interface for TLE validation.

The Stage-5 estimate must remain **independent of the TLE trajectory**. TLEs must not be used as an initial-state prior, position prior, velocity prior, process model or regularisation target.

---

# 2. Scientific Motivation

Stage-4 positions are independent epoch-by-epoch estimates. Their uncertainty is generally anisotropic, especially along the mean line of sight. Stage 5 must exploit temporal continuity and orbital dynamics rather than treating those points as independent.

```text
Stage-4 positions + full covariance
           |
           v
state-space model
           |
           v
forward filter
           |
           v
backward smoother
           |
           v
TLE-independent trajectory estimate
```

The Stage-4 covariance must enter explicitly into the measurement likelihood / Kalman update.

---

# 3. Stage Boundary

```text
Stage 1: common visibility prediction
Stage 2: observation retrieval / validation
Stage 3: station direction estimation
Stage 4: independent 3-D triangulation
Stage 5: statistical trajectory tracking      <-- IMPLEMENT HERE
Stage 6: TLE comparison and residual analysis
Stage 7: station/systematic/anomaly analysis
```

Stage 5 must not implement TLE residuals, TLE correction, station-bias estimation or anomaly attribution.

---

# 4. Stage-4 Inputs

Primary inputs:

```text
positions.parquet
stage4_manifest.json
```

Optional diagnostic input:

```text
station_contributions.parquet
```

Required Stage-4 fields:

```text
time_utc
time_mjd
norad_id
x_ecef_m
y_ecef_m
z_ecef_m
cov_xx
cov_xy
cov_xz
cov_yy
cov_yz
cov_zz
sigma_x_m
sigma_y_m
sigma_z_m
n_stations
condition_number
min_crossing_angle_deg
max_crossing_angle_deg
chi2
reduced_chi2
fit_success
quality_flags
```

Reconstruct the full measurement covariance

\[
R_k =
\begin{bmatrix}
\mathrm{cov}_{xx} & \mathrm{cov}_{xy} & \mathrm{cov}_{xz} \\
\mathrm{cov}_{xy} & \mathrm{cov}_{yy} & \mathrm{cov}_{yz} \\
\mathrm{cov}_{xz} & \mathrm{cov}_{yz} & \mathrm{cov}_{zz}
\end{bmatrix}.
\]

Do not reduce the covariance to diagonal terms.

---

# 5. Coordinate Frames

Stage 4 uses ITRS/ECEF. For the dynamical model, Stage 5 should propagate in an approximately inertial frame.

Use:

```text
GCRS / inertial
```

Processing flow:

```text
Stage-4 ECEF position + covariance
        |
        v
ECEF -> GCRS transformation at each epoch
        |
        v
filter / smoother in inertial coordinates
        |
        v
GCRS state trajectory
        |
        +--> transform position back to ECEF for output
```

Document the exact Astropy/Skyfield frame implementation.

Do not propagate a Newtonian two-body model directly in rotating ECEF coordinates unless all rotating-frame forces are explicitly modelled.

---

# 6. Time Handling

Use actual Stage-4 timestamps. Do not assume uniform cadence.

For adjacent measurements:

\[
\Delta t_k=t_{k+1}-t_k.
\]

Reject or flag:

```text
non-monotonic timestamps
duplicate timestamps
large unexplained gaps
```

Support:

```text
--max-gap-s
```

A long gap must increase prediction uncertainty appropriately.

---

# 7. State Vector

Primary state:

\[
\mathbf x=[x,y,z,v_x,v_y,v_z]^\mathsf T.
\]

Units:

```text
position: m
velocity: m/s
```

Frame:

```text
GCRS / inertial
```

State covariance:

\[
P\in\mathbb R^{6\times6}.
\]

Persist the full covariance.

---

# 8. Dynamical Models

Implement at least two models.

## 8.1 Constant-velocity baseline

\[
\mathbf r_{k+1}=\mathbf r_k+\mathbf v_k\Delta t
\]

\[
\mathbf v_{k+1}=\mathbf v_k.
\]

CLI:

```text
--dynamics constant-velocity
```

This is required as a diagnostic baseline.

## 8.2 Two-body Earth gravity model

Primary production model:

\[
\dot{\mathbf r}=\mathbf v
\]

\[
\dot{\mathbf v}=-\mu\frac{\mathbf r}{|\mathbf r|^3}.
\]

CLI:

```text
--dynamics two-body
```

Default:

```text
two-body
```

Use a documented Earth gravitational parameter from a verified Astropy/Skyfield/geodesy constant source. Avoid unexplained magic constants.

Use an ODE integrator such as `scipy.integrate.solve_ivp` with explicit tolerances.

---

# 9. Optional J2 Model

A `J2` model may be added if cleanly implemented and tested:

```text
--dynamics j2
```

Do not make J2 a prerequisite for initial Stage-5 completion. Two-body dynamics are the initial production model for the short TART observation windows.

---

# 10. Bayesian State-Space Model

State transition:

\[
\mathbf x_{k+1}=f(\mathbf x_k,\Delta t_k)+\mathbf w_k,
\]

with

\[
\mathbf w_k\sim\mathcal N(0,Q_k).
\]

Measurement:

\[
\mathbf y_k=h(\mathbf x_k)+\mathbf v_k,
\]

with

\[
\mathbf v_k\sim\mathcal N(0,R_k).
\]

Preferred implementation: transform Stage-4 ECEF positions and covariances into GCRS before filtering so that

\[
h(\mathbf x)=
\begin{bmatrix}
I_3 & 0
\end{bmatrix}\mathbf x.
\]

---

# 11. Measurement Covariance Transformation

Transform the Stage-4 covariance consistently with the coordinate transformation.

For local Jacobian `J`:

\[
R_{\mathrm{GCRS}}=J R_{\mathrm{ECEF}}J^\mathsf T.
\]

If the position transform is a pure rotation at fixed epoch, use the corresponding rotation matrix.

Validate numerically. Do not rotate positions while leaving covariance in the old frame.

---

# 12. Extended Kalman Filter

Implement a forward EKF.

Prediction:

\[
P_{k|k-1}=F_kP_{k-1|k-1}F_k^\mathsf T+Q_k.
\]

Innovation:

\[
\boldsymbol\nu_k=\mathbf y_k-H_k\hat{\mathbf x}_{k|k-1}.
\]

Innovation covariance:

\[
S_k=H_kP_{k|k-1}H_k^\mathsf T+R_k.
\]

Kalman gain:

\[
K_k=P_{k|k-1}H_k^\mathsf TS_k^{-1}.
\]

Use stable linear solves rather than explicit matrix inversion.

For nonlinear dynamics, propagate or numerically estimate the state-transition Jacobian.

---

# 13. Rauch-Tung-Striebel Smoother

Implement a backward RTS smoother after the forward filter.

Produce both:

```text
filtered state
smoothed state
```

The smoothed trajectory is the primary scientific Stage-5 product for offline analysis.

---

# 14. Process Noise

Support an explicit acceleration-noise parameter:

```text
--process-accel-sigma
```

Units:

```text
m/s^2
```

For constant velocity, use the standard continuous white-acceleration covariance:

\[
Q(\Delta t)=\sigma_a^2
\begin{bmatrix}
\frac{\Delta t^4}{4}I_3 & \frac{\Delta t^3}{2}I_3 \\
\frac{\Delta t^3}{2}I_3 & \Delta t^2 I_3
\end{bmatrix}.
\]

For orbital dynamics, process noise represents unmodelled acceleration / dynamics mismatch.

Record the exact value in the Stage-5 manifest.

Do not tune process noise against TLE agreement.

---

# 15. TLE-Independent Initial State

Initialisation must use Stage-4 data only.

Position:

```text
first valid Stage-4 measurement
```

Velocity hierarchy:

1. weighted local linear fit to the first several Stage-4 positions;
2. finite-difference fallback;
3. zero/minimum-norm velocity only as a flagged last resort.

Support:

```text
--initial-velocity-window 30
```

seconds.

Propagate initial velocity uncertainty into the initial 6x6 covariance.

---

# 16. Innovation Gating

For each measurement calculate

\[
D_k^2=\boldsymbol\nu_k^\mathsf TS_k^{-1}\boldsymbol\nu_k.
\]

For a 3-D Gaussian position measurement:

\[
D_k^2\sim\chi^2_3.
\]

Support:

```text
--innovation-gate-probability
```

Default:

```text
0.9973
```

A gated measurement must remain in diagnostic output with:

```text
INNOVATION_OUTLIER
```

Do not silently discard it.

---

# 17. Robust Measurement Option

Support:

```text
--robust-update none
--robust-update inflate
```

For `inflate`, an outlying measurement may receive inflated covariance rather than full rejection.

Default:

```text
none
```

Keep the first Stage-5 implementation statistically transparent.

---

# 18. Stage-4 Quality Policy

Hard reject by default:

```text
TRIANGULATION_FAILED
INSUFFICIENT_STATIONS
NEGATIVE_RANGE
COVARIANCE_SINGULAR
```

Advisory/configurable:

```text
TIME_MISMATCH
POOR_GEOMETRY
HIGH_LINE_RESIDUAL
HIGH_CHI2
INPUT_DIRECTION_WARNING
CALIBRATION_WARNING
STAGE3_QUALITY_WARNING
BOOTSTRAP_UNSTABLE
```

Support:

```text
--quality-policy strict
--quality-policy permissive
```

Default:

```text
strict
```

---

# 19. Covariance Integrity

Validate every Stage-4 covariance before filtering:

```text
finite
symmetric
positive semidefinite / positive definite within tolerance
reasonable condition number
```

Small numerical asymmetry may be corrected via

\[
R\leftarrow\frac{R+R^\mathsf T}{2}.
\]

Do not silently repair physically invalid covariance matrices.

---

# 20. Numerically Stable Covariance Update

Prefer the Joseph-form covariance update:

\[
P_{k|k}=(I-KH)P_{k|k-1}(I-KH)^\mathsf T+KRK^\mathsf T.
\]

After every prediction/update:

```text
symmetrise covariance
check finite values
check eigenvalues
```

Do not allow silent covariance collapse.

---

# 21. Innovation Statistics

Store for each Stage-4 measurement:

```text
innovation_x_m
innovation_y_m
innovation_z_m
nis
innovation_gate_threshold
innovation_accepted
```

Also compute:

```text
mean NIS
median NIS
fraction above gate
innovation autocorrelation
```

These are TLE-independent model-consistency diagnostics.

---

# 22. Velocity Outputs

Primary velocity outputs:

```text
vx_gcrs_m_s
vy_gcrs_m_s
vz_gcrs_m_s
speed_m_s
```

Do not compare these with TLE velocity until Stage 6.

---

# 23. Acceleration Diagnostic

Record the model acceleration magnitude:

```text
acceleration_m_s2
```

For two-body dynamics:

\[
|\mathbf a|=\mu/r^2.
\]

This is diagnostic only.

---

# 24. Stage-5 Output Layout

```text
results/stage5/
└── YYYYMMDD/
    └── NORAD_<id>/
        ├── trajectory.parquet
        ├── trajectory.csv
        ├── innovations.parquet
        ├── stage5_manifest.json
        ├── diagnostics/
        │   ├── position_track.png
        │   ├── velocity_track.png
        │   ├── uncertainty_track.png
        │   ├── innovations.png
        │   ├── nis.png
        │   ├── stage4_vs_smoothed.png
        │   └── covariance_eigenvalues.png
        └── logs/
            └── stage5.log
```

---

# 25. Trajectory Output Schema

At minimum:

```text
time_utc
time_mjd
norad_id

x_gcrs_filtered_m
y_gcrs_filtered_m
z_gcrs_filtered_m
vx_gcrs_filtered_m_s
vy_gcrs_filtered_m_s
vz_gcrs_filtered_m_s

x_gcrs_smoothed_m
y_gcrs_smoothed_m
z_gcrs_smoothed_m
vx_gcrs_smoothed_m_s
vy_gcrs_smoothed_m_s
vz_gcrs_smoothed_m_s

x_ecef_smoothed_m
y_ecef_smoothed_m
z_ecef_smoothed_m
speed_smoothed_m_s

sigma_x_smoothed_m
sigma_y_smoothed_m
sigma_z_smoothed_m
sigma_vx_smoothed_m_s
sigma_vy_smoothed_m_s
sigma_vz_smoothed_m_s

stage4_x_ecef_m
stage4_y_ecef_m
stage4_z_ecef_m
measurement_used
quality_flags
```

Persist the full filtered and smoothed 6x6 covariance.

---

# 26. Innovation Output

Create:

```text
innovations.parquet
```

Fields:

```text
time_utc
innovation_x_m
innovation_y_m
innovation_z_m
innovation_cov_xx
innovation_cov_xy
innovation_cov_xz
innovation_cov_yy
innovation_cov_yz
innovation_cov_zz
nis
gate_threshold
accepted
stage4_quality_flags
stage5_quality_flags
```

---

# 27. Stage-5 Quality Flags

Suggested flags:

```text
GOOD
INVALID_STAGE4_COVARIANCE
INITIAL_VELOCITY_WEAK
LARGE_TIME_GAP
INNOVATION_OUTLIER
HIGH_NIS
DYNAMICS_MISMATCH
FILTER_DIVERGENCE
COVARIANCE_INVALID
SMOOTHER_FAILED
INSUFFICIENT_MEASUREMENTS
STAGE4_QUALITY_WARNING
```

Support multiple simultaneous flags.

---

# 28. Stage-5 Manifest

Write:

```text
stage5_manifest.json
```

with:

```json
{
  "schema_version": 1,
  "stage": 5,
  "input_stage4_manifest": "...",
  "input_sha256": "...",
  "norad_id": "...",
  "state_frame": "GCRS",
  "measurement_frame": "ITRS/ECEF",
  "state_model": {},
  "filter": {},
  "smoother": {},
  "process_noise": {},
  "innovation_gate": {},
  "quality_policy": {},
  "software": {},
  "provenance": {},
  "outputs": {}
}
```

Record exact inputs, git commit, command line, dependency versions, dynamics, integrator tolerances, process noise, initialisation method, innovation gate and quality policy.

---

# 29. CLI

Implement:

```bash
tart-track
```

Example:

```bash
tart-track \
    --stage4-manifest results/stage4/20260625/NORAD_40534/stage4_manifest.json \
    --dynamics two-body \
    --output-dir results/stage5
```

Support:

```text
--stage4-manifest
--output-dir
--dynamics
--process-accel-sigma
--initial-velocity-window
--innovation-gate-probability
--quality-policy
--max-gap-s
--robust-update
--no-smoother
--no-plots
--verbose
```

Default dynamics:

```text
two-body
```

---

# 30. Repository Integration

Add:

```text
src/tart_gnss/stage5/
├── __init__.py
├── cli.py
├── loader.py
├── frames.py
├── dynamics.py
├── initialisation.py
├── filter.py
├── smoother.py
├── innovations.py
├── quality.py
├── output.py
├── plotting.py
└── errors.py
```

Tests:

```text
tests/stage5/
├── __init__.py
├── test_loader.py
├── test_frames.py
├── test_dynamics.py
├── test_initialisation.py
├── test_filter.py
├── test_smoother.py
├── test_innovations.py
├── test_quality.py
├── test_output.py
├── test_cli.py
├── test_synthetic_tracking.py
└── test_covariance_consistency.py
```

Documentation:

```text
docs/stage5_tracking.md
```

---

# 31. API Design

Suggested interfaces:

```python
def load_stage4_positions(manifest_path: Path) -> Stage4TrajectoryInput:
    """Load validated Stage-4 positions, covariances and provenance."""
```

```python
def ecef_measurement_to_gcrs(
    position_ecef_m: np.ndarray,
    covariance_ecef_m2: np.ndarray,
    obstime: datetime,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform one ECEF position and covariance into the inertial frame."""
```

```python
def propagate_two_body(
    state: np.ndarray,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate position/velocity and return the state-transition Jacobian."""
```

```python
def ekf_filter(
    measurements: list[TrajectoryMeasurement],
    config: FilterConfig,
) -> FilterResult:
    """Run forward extended Kalman filtering over Stage-4 measurements."""
```

```python
def rts_smooth(filtered: FilterResult) -> SmoothedResult:
    """Run backward Rauch-Tung-Striebel smoothing."""
```

---

# 32. Synthetic Validation

Required tests:

## 32.1 Exact two-body trajectory

Generate a GNSS-like orbit segment and noisy Stage-4-like positions. Recover position and velocity.

## 32.2 Strong covariance anisotropy

Use small transverse uncertainty and large LOS uncertainty. Verify the full covariance changes the solution correctly.

## 32.3 Irregular cadence

Perturb observation times. Verify actual `dt` is used.

## 32.4 Missing epochs

Remove measurements. Verify uncertainty grows during prediction-only intervals.

## 32.5 Outlier

Inject one bad Stage-4 measurement. Verify innovation gating detects it.

## 32.6 Process-noise mismatch

Generate small unmodelled acceleration. Verify process noise affects innovation consistency as expected.

## 32.7 No-TLE test

The entire Stage-5 synthetic pipeline must run with no TLE data available. This is mandatory.

---

# 33. Monte Carlo Coverage Validation

Across repeated synthetic realisations, evaluate state coverage.

At minimum:

```text
fraction of true positions within 1-sigma ellipsoid
fraction within 2-sigma ellipsoid
fraction within 3-sigma ellipsoid
velocity coverage
```

Large miscalibration must be documented and fixed or explicitly flagged.

---

# 34. NEES Validation

For synthetic truth only:

\[
\mathrm{NEES}=
(\hat{\mathbf x}-\mathbf x_{true})^\mathsf T
P^{-1}
(\hat{\mathbf x}-\mathbf x_{true}).
\]

For a six-dimensional Gaussian state:

\[
\mathrm{NEES}\sim\chi^2_6.
\]

Do not compute NEES on real data by treating TLE positions as truth.

---

# 35. Dynamics-Model Comparison

Compare:

```text
constant-velocity
two-body
```

using only TLE-independent diagnostics:

```text
innovation RMS
median NIS
innovation autocorrelation
accepted measurement fraction
```

Do not select a dynamics model based on TLE agreement.

---

# 36. Real-Data Acceptance Test

Use the existing Stage-4 result for:

```text
2026-06-25
NORAD 40534
```

if available locally.

Run:

```bash
tart-track \
    --stage4-manifest results/stage4/20260625/NORAD_40534/stage4_manifest.json \
    --dynamics two-body \
    --output-dir results/stage5
```

Inspect:

```text
trajectory continuity
filtered vs smoothed position
velocity continuity
state uncertainty
NIS over time
innovation autocorrelation
Stage-4 residuals around smoothed trajectory
outlier epochs
gap behaviour
```

Do not interpret TLE residuals yet.

---

# 37. Diagnostic Plots

Implement at minimum:

1. Stage-4 measurements vs filtered and smoothed `x(t), y(t), z(t)`;
2. `vx(t), vy(t), vz(t), speed(t)`;
3. innovation components;
4. NIS plus gate threshold;
5. position uncertainty vs time;
6. covariance eigenvalues vs time.

Use matplotlib.

---

# 38. Logging

At INFO level report:

```text
input Stage-4 manifest
number of input epochs
number retained after quality filtering
time span
median cadence
maximum gap
dynamics model
process acceleration sigma
initial position
initial velocity
innovation gate
number of rejected innovations
median NIS
filter success
smoother success
median filtered position uncertainty
median smoothed position uncertainty
output paths
```

Use `logging`.

---

# 39. Failure Behaviour

Add explicit exceptions:

```text
Stage5Error
InvalidStage4ManifestError
MissingPositionDataError
InvalidCovarianceError
FrameTransformError
InitialisationError
DynamicsError
FilterDivergenceError
SmootherError
```

CLI exit codes:

```text
0 successful scientific run
2 invalid input/configuration
3 completed but insufficient valid trajectory estimates
1 unexpected execution failure
```

---

# 40. Scientific Guardrails

Stage 5 must not create artificial precision.

Mandatory rules:

1. retain full Stage-4 covariance;
2. make process noise explicit;
3. make gaps increase uncertainty;
4. preserve rejected measurements in diagnostics;
5. report smoothed covariance;
6. keep TLEs out of estimation;
7. never rescale covariance merely to improve later TLE agreement.

---

# 41. Acceptance Criteria

## CLI

- [ ] `tart-track` exists.
- [ ] Existing `tart-*` commands remain unchanged.
- [ ] `tests/test_cli_names.py` includes `tart-track`.
- [ ] No non-`tart-*` Stage-5 public entry point exists.

## Functional

- [ ] Stage-4 positions and full covariance are loaded.
- [ ] ECEF positions/covariances transform to inertial coordinates.
- [ ] Irregular timestamps are supported.
- [ ] Constant-velocity dynamics work.
- [ ] Two-body dynamics work.
- [ ] EKF filtering works.
- [ ] RTS smoothing works.
- [ ] Innovation gating works.
- [ ] Full 6x6 covariance is persisted.
- [ ] Innovations are persisted.
- [ ] Stage-5 manifest is written.
- [ ] Diagnostic plots are produced.

## Scientific

- [ ] Synthetic trajectory is recovered accurately.
- [ ] Missing measurements increase uncertainty.
- [ ] Outliers are detected statistically.
- [ ] Full anisotropic Stage-4 covariance affects the estimate correctly.
- [ ] Monte Carlo coverage is approximately calibrated.
- [ ] Synthetic NEES is reasonable.
- [ ] The pipeline runs without TLE data.
- [ ] Two-body dynamics outperform constant velocity on orbital synthetic data using TLE-independent metrics.

## Engineering

- [ ] Stage-1 tests pass.
- [ ] Stage-2 tests pass.
- [ ] Stage-3 tests pass.
- [ ] Stage-4 tests pass.
- [ ] Stage-5 tests pass.
- [ ] Registry tests pass.
- [ ] No network required for unit tests.
- [ ] Large real observational data are not committed.
- [ ] Public functions have docstrings.
- [ ] Logging is consistent.

---

# 42. Do Not Implement Yet

Do not implement in Stage 5:

```text
TLE residuals
TLE bias correction
radial/along-track/cross-track TLE errors
orbit anomaly classification
station bias estimation
calibration correction
RFI attribution
automated sensor-health scoring
```

These belong to Stages 6 and 7.

---

# 43. Expected End State

Public pipeline:

```bash
tart-stations-update ...
tart-stations-list ...
tart-common-gnss ...
tart-retrieve ...
tart-localise ...
tart-triangulate ...
tart-track ...
```

Conceptually:

```text
TART visibilities
      |
tart-localise
      |
station directions + covariance
      |
tart-triangulate
      |
independent 3-D positions + covariance
      |
tart-track
      |
smoothed TLE-independent trajectory
position + velocity + full covariance
      |
Stage 6
      |
comparison with TLE prediction
```

---

# 44. Stage-6 Interface

Stage 5 must expose, at each epoch:

```text
time
smoothed GCRS position
smoothed GCRS velocity
smoothed ECEF position
full 6x6 state covariance
innovation/measurement quality
```

Stage 6 will independently propagate the archived TLE to these epochs and compute observed-minus-predicted residuals.

Do not precompute TLE residuals in Stage 5.

---

# 45. Final Deliverables

The coding agent must deliver:

```text
src/tart_gnss/stage5/*
tests/stage5/*
docs/stage5_tracking.md
README.md update
pyproject.toml entry-point/dependency updates
tests/test_cli_names.py update
tart_gnss_stage5_implementation.md
```

Also provide:

```text
synthetic Stage-5 example
real-data Stage-5 example if local Stage-4 data are available
```

Before completion run:

```bash
pytest
```

and report:

```text
passed
failed
skipped
```

Run an end-to-end synthetic demonstration and report:

```text
true initial position
true initial velocity
recovered smoothed position
recovered smoothed velocity
position RMSE
velocity RMSE
mean / median NIS
innovation rejection fraction
predicted state uncertainty
empirical Monte Carlo error
constant-velocity diagnostics
two-body diagnostics
```

The implementation is complete only when the Stage-5 trajectory is demonstrably derived from Stage-4 measurements and their uncertainties, with no dependence on TLE predictions.
