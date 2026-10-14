# TART-GNSS Stage 6 Implementation Specification
## TLE Trajectory Validation and Residual Analysis

**Repository:** `git@github.com:ifornax/tart-gnss.git`  
**Base branch:** `main`  
**Stage status before this work:** Stages 1–5 implemented  
**This specification implements:** Stage 6 only  
**Canonical Stage-6 CLI:** `tart-validate`

---

# 0. Permanent CLI Naming Convention

Every public command in this repository must begin with `tart-`.

Current canonical commands:

```text
tart-stations-update
tart-stations-list
tart-common-gnss
tart-retrieve
tart-localise
tart-triangulate
tart-track
```

Stage 6 must add:

```text
tart-validate
```

Do not introduce alternative public command names. Update `pyproject.toml`, `README.md`, documentation, examples, and `tests/test_cli_names.py`. Historical provenance files must remain unchanged.

---

# 1. Objective

Implement **Stage 6: comparison of the TLE-independent measured trajectory from Stage 5 against a TLE-predicted satellite trajectory**.

Stage 5 provides a smoothed state

\[
\hat{\mathbf x}_{\rm obs}(t)=
\begin{bmatrix}
\hat{\mathbf r}_{\rm obs}(t)\\
\hat{\mathbf v}_{\rm obs}(t)
\end{bmatrix}
\]

with full covariance

\[
P_{\rm obs}(t)\in\mathbb R^{6\times6}.
\]

Stage 6 must independently propagate the selected TLE to the exact Stage-5 epochs:

\[
\mathbf x_{\rm TLE}(t)=
\begin{bmatrix}
\mathbf r_{\rm TLE}(t)\\
\mathbf v_{\rm TLE}(t)
\end{bmatrix}.
\]

Then calculate

\[
\Delta\mathbf r(t)=\mathbf r_{\rm obs}(t)-\mathbf r_{\rm TLE}(t)
\]

and

\[
\Delta\mathbf v(t)=\mathbf v_{\rm obs}(t)-\mathbf v_{\rm TLE}(t).
\]

**Permanent residual sign convention:**

```text
observed - TLE
```

This must be documented and unit-tested.

Stage 6 is the first stage where the TLE enters the scientific comparison.

---

# 2. Scientific Questions

Stage 6 must determine:

1. the 3-D separation between the measured and TLE-predicted trajectories;
2. whether the mismatch is primarily radial, along-track, or cross-track;
3. how residuals evolve with time;
4. whether residuals are large relative to the **measurement uncertainty** from Stage 5;
5. whether the selected TLE is sufficiently contemporaneous for precision validation;
6. whether nearby element sets imply significant reference-model dispersion;
7. whether the residuals show a constant bias, trend, or time-correlated structure.

Do not attribute residuals to a physical anomaly in Stage 6. Attribution belongs to Stage 7.

---

# 3. Critical Statistical Principle: TLE Is Not Exact Ground Truth

A TLE normally provides no formal state covariance. Stage 6 must therefore never report a formal “sigma discrepancy” while silently treating the TLE as exact.

## 3.1 Measurement-conditional statistic

Using only Stage-5 position covariance:

\[
D_{\rm obs}^2=
\Delta\mathbf r^\mathsf T C_{\rm obs}^{-1}\Delta\mathbf r.
\]

Output:

```text
position_mahalanobis_obs_cov
```

Meaning: discrepancy relative to TART measurement uncertainty **conditional on the selected TLE being fixed**.

For the full 6-D state:

```text
state_mahalanobis_obs_cov
```

## 3.2 Reference-aware statistic

If an empirical TLE-reference covariance is available from several nearby element sets:

\[
C_{\rm total}=C_{\rm obs}+C_{\rm TLE,empirical}.
\]

Then calculate:

```text
position_mahalanobis_total_cov
state_mahalanobis_total_cov
```

The TLE covariance must be labelled:

```text
empirical reference dispersion
```

and never called a formal orbit-determination covariance.

---

# 4. Stage-5 Inputs

Primary inputs:

```text
stage5_manifest.json
trajectory.parquet
```

Optional:

```text
innovations.parquet
```

Required trajectory fields include:

```text
time_utc
time_mjd
norad_id
x_gcrs_smoothed_m
y_gcrs_smoothed_m
z_gcrs_smoothed_m
vx_gcrs_smoothed_m_s
vy_gcrs_smoothed_m_s
vz_gcrs_smoothed_m_s
x_ecef_smoothed_m
y_ecef_smoothed_m
z_ecef_smoothed_m
valid
quality_flags
```

Stage 6 must reconstruct the full smoothed 6×6 covariance from the Stage-5 covariance columns. Do not reduce it to diagonal terms.

---

# 5. TLE Source Hierarchy

Support three reference modes.

## 5.1 Exact archived Stage-1 TLE

Use the exact archived Stage-1 element set when scientifically appropriate.

CLI:

```text
--tle-source stage1
```

## 5.2 Explicit TLE file

Preferred for historical precision validation:

```bash
--tle-file path/to/catalogue.tle
```

## 5.3 Explicit download

If live download is supported, it must be explicitly requested and the downloaded file must be archived with provenance.

Never silently download a current TLE and use it to validate historical observations.

---

# 6. TLE Epoch and Age

For every Stage-5 epoch calculate

\[
\Delta t_{\rm TLE}=t_{\rm obs}-t_{\rm TLEepoch}.
\]

Output:

```text
tle_epoch_utc
tle_age_days
tle_propagation_direction
```

Convention:

```text
positive age = propagation after TLE epoch
negative age = propagation before TLE epoch
```

Flags:

```text
TLE_AGE_WARNING
TLE_STALE
TLE_BACKPROPAGATED
```

Suggested defaults:

```text
--tle-age-warning-days 3
--tle-age-stale-days 14
```

These are pipeline quality thresholds, not universal SGP4 accuracy limits.

If the TLE is stale, still allow descriptive outputs but set:

```text
precision_validation = false
```

The existing 2026-06-25 example with a TLE propagated about 103 days backward must therefore be labelled unsuitable for precision validation.

---

# 7. TLE Propagation

Use SGP4 through the project's existing Skyfield/sgp4 infrastructure where practical.

Propagate the selected TLE at the **exact Stage-5 epochs**.

Do not interpolate a sparse TLE trajectory when direct propagation is available.

Output the predicted state in:

```text
GCRS-compatible inertial coordinates
ITRS/ECEF coordinates
```

Units:

```text
metres
metres/second
```

---

# 8. Frame Consistency

Primary comparison frame:

```text
GCRS
```

Secondary convenience frame:

```text
ITRS/ECEF position
```

Observed and predicted states must be compared at the same epoch and in the same frame.

Never compare TEME directly with GCRS or ECEF velocity directly with inertial velocity without the correct transformation.

Add numerical frame-consistency tests.

---

# 9. TEME / SGP4 Handling

SGP4 naturally produces a TEME state. Document explicitly how that state is transformed into the Stage-5 comparison frame.

Documentation must state:

```text
SGP4 source frame
transformation library
comparison frame
time scale
Earth-orientation assumptions
```

Do not assume TEME equals GCRS/J2000.

---

# 10. Cartesian Residuals

Calculate in GCRS:

```text
dr_x_m
dr_y_m
dr_z_m
position_residual_3d_m

dv_x_m_s
dv_y_m_s
dv_z_m_s
velocity_residual_3d_m_s
```

where

\[
\|\Delta\mathbf r\|
\]

and

\[
\|\Delta\mathbf v\|
\]

are the 3-D position and velocity residual magnitudes.

---

# 11. RTN / RIC Orbital Frame

Use the TLE-predicted state to define the reference RTN frame.

Radial:

\[
\hat{\mathbf R}=\frac{\mathbf r_{\rm TLE}}{\|\mathbf r_{\rm TLE}\|}.
\]

Cross-track / normal:

\[
\hat{\mathbf N}=
\frac{\mathbf r_{\rm TLE}\times\mathbf v_{\rm TLE}}
{\|\mathbf r_{\rm TLE}\times\mathbf v_{\rm TLE}\|}.
\]

Along-track / transverse:

\[
\hat{\mathbf T}=\hat{\mathbf N}\times\hat{\mathbf R}.
\]

Position residuals:

```text
dr_radial_m
dr_along_track_m
dr_cross_track_m
```

Velocity residuals:

```text
dv_radial_m_s
dv_along_track_m_s
dv_cross_track_m_s
```

Document:

```text
RTN = RIC convention used here
R = radial
T/I = along-track / in-track
N/C = cross-track
```

---

# 12. Covariance Rotation into RTN

Let the Stage-5 smoothed GCRS position covariance be `C_r`.

Define

\[
Q=
\begin{bmatrix}
\hat R^\mathsf T\\
\hat T^\mathsf T\\
\hat N^\mathsf T
\end{bmatrix}.
\]

Then

\[
C_{RTN}=Q C_r Q^\mathsf T.
\]

Output:

```text
sigma_radial_obs_m
sigma_along_track_obs_m
sigma_cross_track_obs_m
cov_radial_along_obs_m2
cov_radial_cross_obs_m2
cov_along_cross_obs_m2
```

Use the full covariance, not only Cartesian standard deviations.

---

# 13. Velocity Covariance

Extract the velocity block from the full Stage-5 smoothed 6×6 covariance and rotate it into RTN.

Output:

```text
sigma_dv_radial_obs_m_s
sigma_dv_along_track_obs_m_s
sigma_dv_cross_track_obs_m_s
```

Retain the full state covariance for advanced analyses.

---

# 14. Timing-Like Along-Track Diagnostic

Calculate:

\[
\Delta t_{\rm along}=\frac{\Delta r_T}{\|\mathbf v_{\rm TLE}\|}.
\]

Output:

```text
equivalent_along_track_time_offset_s
```

This is a diagnostic only. Do not interpret it automatically as a clock error.

---

# 15. Angular and Range Residuals

Calculate geocentric angular separation:

\[
\Delta\theta=\cos^{-1}(\hat{\mathbf r}_{\rm obs}\cdot\hat{\mathbf r}_{\rm TLE}).
\]

Output:

```text
geocentric_angular_residual_deg
```

Calculate geocentric range residual:

\[
\Delta\rho=\|\mathbf r_{\rm obs}\|-\|\mathbf r_{\rm TLE}\|.
\]

Output:

```text
geocentric_range_residual_m
```

---

# 16. TLE Ensemble / Empirical Reference Dispersion

Support optional reference ensembles:

```text
--tle-ensemble-dir
--tle-ensemble-files file1.tle file2.tle ...
--tle-ensemble-window-days
```

Suggested default ensemble window:

```text
7 days
```

For each TLE of the same NORAD ID:

1. propagate to each Stage-5 epoch;
2. transform into the same frame;
3. compute the propagated ensemble mean;
4. calculate empirical state covariance.

Output:

```text
n_tle_ensemble
tle_position_dispersion_covariance
tle_velocity_dispersion_covariance
```

This is empirical reference dispersion, not formal SGP4 covariance.

The explicitly selected TLE remains the primary reference. Do not silently replace it with the ensemble mean.

---

# 17. Measurement-Conditional Mahalanobis Distances

Position:

\[
D_{\rm obs,pos}^2=
\Delta\mathbf r^\mathsf T C_r^{-1}\Delta\mathbf r.
\]

State:

\[
D_{\rm obs,state}^2=
\Delta\mathbf x^\mathsf T P_{\rm obs}^{-1}\Delta\mathbf x.
\]

Output:

```text
position_mahalanobis_obs_cov
state_mahalanobis_obs_cov
```

These are conditional on the selected TLE state being treated as fixed.

---

# 18. Reference-Aware Mahalanobis Distances

When empirical TLE covariance is available:

\[
C_{\rm tot}=C_{\rm obs}+C_{\rm TLE,emp}.
\]

Output:

```text
position_mahalanobis_total_cov
state_mahalanobis_total_cov
```

Require at least three ensemble TLEs before producing the total-covariance statistic.

Flag small ensembles:

```text
TLE_ENSEMBLE_SMALL
```

---

# 19. Time-Series Summary Statistics

For RTN position residuals report:

```text
mean
median
RMS
MAD-based robust scale
16th/84th percentiles
minimum
maximum
```

for:

```text
radial
along-track
cross-track
```

Do not assume residual epochs are statistically independent.

---

# 20. Temporal Trend Analysis

Fit descriptive linear trends:

\[
\Delta r_j(t)=a_j+b_j(t-t_0)
\]

for radial, along-track and cross-track components.

Use weighted least squares with projected Stage-5 covariance where practical.

Output:

```text
residual_slope_radial_m_s
residual_slope_along_track_m_s
residual_slope_cross_track_m_s
```

These are mismatch rates, not orbit corrections.

---

# 21. Residual Autocorrelation

Calculate at minimum lag-1 autocorrelation of:

```text
radial residual
along-track residual
cross-track residual
```

Strong autocorrelation indicates structured residuals rather than independent measurement noise.

Do not attribute the cause in Stage 6.

---

# 22. Stage-5 Quality Policy

Default strict policy excludes rows with blocking Stage-5 conditions:

```text
NOT_INITIALISED
SMOOTHER_FAILED
COVARIANCE_INVALID
INSUFFICIENT_MEASUREMENTS
```

Advisory flags may propagate:

```text
FILTER_DIVERGENCE
INNOVATION_OUTLIER
LARGE_TIME_GAP
EXTRAPOLATED
HIGH_NIS
DYNAMICS_MISMATCH
STAGE4_QUALITY_WARNING
INITIAL_VELOCITY_WEAK
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

# 23. Stage-6 Quality Flags

Suggested flags:

```text
GOOD
INVALID_STAGE5_STATE
INVALID_STAGE5_COVARIANCE
TLE_NOT_FOUND
TLE_NORAD_MISMATCH
TLE_AGE_WARNING
TLE_STALE
TLE_BACKPROPAGATED
TLE_FRAME_FAILURE
TLE_ENSEMBLE_SMALL
TLE_REFERENCE_DISPERSION_HIGH
REFERENCE_UNCERTAINTY_UNKNOWN
HIGH_POSITION_RESIDUAL
HIGH_VELOCITY_RESIDUAL
HIGH_CONDITIONAL_MAHALANOBIS
STRONG_RESIDUAL_TREND
RESIDUAL_AUTOCORRELATED
STAGE5_QUALITY_WARNING
NOT_PRECISION_VALIDATION
```

Support multiple flags simultaneously.

---

# 24. Precision Validation Status

Every run must state:

```text
precision_validation: true/false
```

Suggested minimum conditions for `true`:

1. TLE exists and NORAD matches;
2. TLE age does not exceed `--tle-age-stale-days`;
3. Stage-5 trajectory is valid;
4. frame transformation succeeds;
5. no critical state-covariance failure.

This Boolean does not mean the TLE is exact. It means the comparison meets the pipeline's minimum conditions for precision analysis.

---

# 25. User-Defined Residual Thresholds

Do not hard-code arbitrary kilometre thresholds as scientific truth.

Allow:

```text
--position-residual-warning-km
--velocity-residual-warning-m-s
```

If omitted, do not set high-residual flags based on arbitrary defaults.

---

# 26. Output Layout

```text
results/stage6/
└── <YYYYMMDD>/
    └── NORAD_<id>/
        ├── residuals.parquet
        ├── residuals.csv
        ├── tle_trajectory.parquet
        ├── summary.json
        ├── stage6_manifest.json
        ├── diagnostics/
        │   ├── rtn_residuals.png
        │   ├── position_residual_norm.png
        │   ├── velocity_residuals.png
        │   ├── angular_residual.png
        │   ├── residual_vs_uncertainty.png
        │   ├── mahalanobis.png
        │   ├── tle_age.png
        │   └── trajectory_overlay.png
        └── logs/
            └── stage6.log
```

---

# 27. Residual Output Schema

At minimum:

```text
time_utc
time_mjd
norad_id
tle_epoch_utc
tle_age_days
obs_x_gcrs_m
obs_y_gcrs_m
obs_z_gcrs_m
obs_vx_gcrs_m_s
obs_vy_gcrs_m_s
obs_vz_gcrs_m_s
tle_x_gcrs_m
tle_y_gcrs_m
tle_z_gcrs_m
tle_vx_gcrs_m_s
tle_vy_gcrs_m_s
tle_vz_gcrs_m_s
dr_x_m
dr_y_m
dr_z_m
position_residual_3d_m
dv_x_m_s
dv_y_m_s
dv_z_m_s
velocity_residual_3d_m_s
dr_radial_m
dr_along_track_m
dr_cross_track_m
dv_radial_m_s
dv_along_track_m_s
dv_cross_track_m_s
sigma_radial_obs_m
sigma_along_track_obs_m
sigma_cross_track_obs_m
geocentric_angular_residual_deg
geocentric_range_residual_m
equivalent_along_track_time_offset_s
position_mahalanobis_obs_cov
state_mahalanobis_obs_cov
position_mahalanobis_total_cov
state_mahalanobis_total_cov
n_tle_ensemble
precision_validation
valid
quality_flags
```

Use NaN for total-covariance metrics when no empirical reference covariance exists.

---

# 28. TLE Trajectory Output

Write:

```text
tle_trajectory.parquet
```

including:

```text
time_utc
norad_id
tle_epoch_utc
tle_age_days
x_gcrs_m
y_gcrs_m
z_gcrs_m
vx_gcrs_m_s
vy_gcrs_m_s
vz_gcrs_m_s
x_ecef_m
y_ecef_m
z_ecef_m
RTN basis vectors
```

This avoids re-propagating the TLE for later Stage-7 analysis.

---

# 29. Summary Output

Write `summary.json` containing:

```text
number of epochs
number valid
time span
selected TLE
TLE epoch
median TLE age
maximum absolute TLE age
precision-validation status
3-D position residual mean/median/RMS/max
velocity residual mean/median/RMS/max
RTN mean/median/RMS/MAD
RTN trend slopes
conditional Mahalanobis summaries
ensemble-TLE count and dispersion summary
residual autocorrelation
flag counts
```

---

# 30. Stage-6 Manifest

Write `stage6_manifest.json` with:

```json
{
  "schema_version": 1,
  "stage": 6,
  "input_stage5_manifest": "...",
  "input_stage5_sha256": "...",
  "norad_id": "...",
  "tle_reference": {},
  "tle_ensemble": {},
  "frames": {},
  "residual_definition": "observed_minus_tle",
  "quality_policy": {},
  "precision_validation": {},
  "software": {},
  "provenance": {},
  "outputs": {}
}
```

Record exact TLE lines, file path, SHA-256, source, epoch, age thresholds, SGP4 implementation/version, frame-transform method, Stage-5 input hash, git commit, command line, and software versions.

---

# 31. CLI

Implement:

```bash
tart-validate
```

Example:

```bash
tart-validate \
    --stage5-manifest results/stage5/20260625/NORAD_40534/stage5_manifest.json \
    --tle-file data/tle/gnss_20260625.tle \
    --output-dir results/stage6
```

Support:

```text
--stage5-manifest
--tle-file
--tle-source
--tle-ensemble-files
--tle-ensemble-dir
--tle-ensemble-window-days
--tle-age-warning-days
--tle-age-stale-days
--quality-policy
--position-residual-warning-km
--velocity-residual-warning-m-s
--allow-stale-reference
--no-plots
--verbose
```

Require exactly one primary TLE source.

---

# 32. Repository Integration

Add:

```text
src/tart_gnss/stage6/
├── __init__.py
├── cli.py
├── loader.py
├── tle_reference.py
├── propagation.py
├── frames.py
├── rtn.py
├── residuals.py
├── reference_uncertainty.py
├── statistics.py
├── quality.py
├── output.py
├── plotting.py
└── errors.py
```

Tests:

```text
tests/stage6/
├── __init__.py
├── test_loader.py
├── test_tle_reference.py
├── test_propagation.py
├── test_frames.py
├── test_rtn.py
├── test_residuals.py
├── test_reference_uncertainty.py
├── test_statistics.py
├── test_quality.py
├── test_output.py
├── test_cli.py
└── test_synthetic_validation.py
```

Documentation:

```text
docs/stage6_validation.md
```

---

# 33. API Design

Suggested interfaces:

```python
def load_stage5_trajectory(manifest_path: Path) -> Stage5Trajectory:
    """Load the Stage-5 smoothed trajectory and full covariance."""
```

```python
def load_tle_reference(norad_id: int, tle_file: Path) -> TLEReference:
    """Load the selected NORAD element set with provenance."""
```

```python
def propagate_tle(tle: TLEReference, times: np.ndarray) -> TLETrajectory:
    """Propagate the TLE with SGP4 at exact Stage-5 epochs."""
```

```python
def rtn_basis(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """Return the orthonormal R/T/N basis from a reference state."""
```

```python
def compute_residuals(observed: StateTrajectory, reference: StateTrajectory) -> ResidualTable:
    """Calculate observed-minus-reference Cartesian and RTN residuals."""
```

```python
def empirical_tle_covariance(references: list[TLEReference], times: np.ndarray) -> ReferenceCovariance:
    """Estimate propagated-state dispersion from an ensemble of TLEs."""
```

---

# 34. Synthetic Validation

Mandatory tests:

## 34.1 Zero-residual

Generate a truth trajectory from a known TLE and use the same propagated state as the observed trajectory. Residuals must be approximately zero.

## 34.2 Known Cartesian offset

Inject a known inertial offset and verify the `observed - TLE` sign convention.

## 34.3 Known RTN offset

Inject:

```text
+1000 m radial
+2000 m along-track
-500 m cross-track
```

and verify recovery.

## 34.4 Velocity offset

Inject known velocity perturbations and verify RTN velocity residuals.

## 34.5 Covariance rotation

Verify GCRS-to-RTN covariance rotation numerically.

## 34.6 TLE age

Test fresh, warning-age, stale, and back-propagated references.

## 34.7 No formal reference covariance

With one TLE only:

```text
position_mahalanobis_obs_cov exists
position_mahalanobis_total_cov is NaN
REFERENCE_UNCERTAINTY_UNKNOWN is set
```

## 34.8 TLE ensemble

Use nearby synthetic element sets and verify empirical reference dispersion and total-covariance statistics.

---

# 35. Frame Validation

Frame tests are mandatory.

Validate:

```text
TEME -> comparison frame
GCRS consistency
ITRS position consistency
velocity units/orientation
```

A frame error can mimic a physical residual and is unacceptable.

---

# 36. RTN Basis Validation

Test orthogonality and handedness:

\[
\hat R\cdot\hat T=0,
\quad
\hat R\cdot\hat N=0,
\quad
\hat T\cdot\hat N=0,
\]

and

\[
\det(Q)=+1.
\]

Verify:

```text
R x T = N
```

within numerical tolerance.

---

# 37. Real-Data Acceptance Test

Use Stage-5 result:

```text
2026-06-25
NORAD 40534
```

only if a suitably contemporaneous TLE is available.

If only the historical TLE propagated about 103 days backward is available:

1. run Stage 6 for software validation;
2. set `precision_validation = false`;
3. flag `TLE_STALE`, `TLE_BACKPROPAGATED`, and `NOT_PRECISION_VALIDATION`;
4. do not interpret the residual magnitude as telescope astrometric accuracy.

For the scientific acceptance result, obtain an element set close to the actual observation epoch.

---

# 38. Diagnostic Plots

Implement at minimum:

```text
RTN residuals vs time
3-D position residual norm vs time
RTN velocity residuals vs time
geocentric angular residual vs time
conditional Mahalanobis vs time
reference-aware Mahalanobis vs time when available
observed/TLE trajectory overlay
TLE age vs time
```

Use matplotlib.

---

# 39. Statistical Reporting

Report residuals in physically meaningful units:

```text
m / km
m/s
degrees
seconds-equivalent along-track
```

Report robust and classical summaries:

```text
median
MAD
percentiles
mean
RMS
```

---

# 40. Interpretation Guardrails

A large residual does not automatically imply:

```text
bad TLE
bad telescope
satellite manoeuvre
clock error
calibration error
RFI contamination
```

A small residual does not prove absolute astrometric accuracy or correct uncertainty calibration.

Stage 7 is responsible for systematic attribution.

---

# 41. Do Not Fit Away the Residual

Stage 6 must not:

```text
shift the TLE trajectory to minimise residuals
fit and subtract a time offset
fit and subtract an RTN translation
rescale Stage-5 covariance to force reduced chi-square near 1
adjust TLE orbital elements
```

A descriptive constant bias may be reported separately, but raw observed-minus-TLE residuals remain the primary result.

---

# 42. Reproducibility

A Stage-6 result must be reproducible from:

```text
Stage-5 manifest
Stage-5 trajectory
exact selected TLE
TLE SHA-256
SGP4 implementation/version
frame implementation
Stage-6 manifest
git commit
command line
```

No live network lookup should be required to reproduce a completed run.

---

# 43. Logging

At INFO level report:

```text
Stage-5 input
NORAD ID
number of valid epochs
selected TLE source
TLE epoch
median/max TLE age
forward/backward propagation
precision-validation status
ensemble usage/count
median and RMS 3-D residual
median RTN residuals
median conditional Mahalanobis
output paths
```

Use `logging`.

---

# 44. Failure Behaviour

Add:

```text
Stage6Error
InvalidStage5ManifestError
MissingTrajectoryError
TLEReferenceError
TLENoradMismatchError
TLEPropagationError
FrameTransformError
InvalidStateCovarianceError
RTNGeometryError
ReferenceCovarianceError
```

CLI exit codes:

```text
0 successful comparison
2 invalid input/configuration
3 comparison completed but not scientifically valid for precision analysis
1 unexpected failure
```

`--allow-stale-reference` may permit exit 0 for descriptive analysis while still flagging stale-reference status.

---

# 45. Acceptance Criteria

## CLI

- [ ] `tart-validate` exists.
- [ ] Existing `tart-*` commands remain unchanged.
- [ ] `tests/test_cli_names.py` includes `tart-validate`.

## Functional

- [ ] Stage-5 trajectory and full covariance load.
- [ ] Exact TLE is selected by NORAD ID.
- [ ] TLE propagates at exact Stage-5 epochs.
- [ ] Frames are consistent.
- [ ] Cartesian position and velocity residuals are calculated.
- [ ] RTN position and velocity residuals are calculated.
- [ ] Stage-5 covariance is rotated into RTN.
- [ ] TLE age is calculated.
- [ ] Precision-validation status is produced.
- [ ] Optional TLE ensemble works.
- [ ] Residual and TLE trajectory tables are written.
- [ ] Summary and manifest are written.
- [ ] Diagnostic plots are generated.

## Scientific

- [ ] Residual sign is `observed - TLE`.
- [ ] Zero-residual synthetic case returns approximately zero.
- [ ] Known RTN injections are recovered.
- [ ] Velocity residuals are correctly oriented.
- [ ] Single-TLE runs never pretend to have formal reference covariance.
- [ ] Conditional Mahalanobis is labelled correctly.
- [ ] Ensemble dispersion is labelled empirical.
- [ ] Stale TLE is clearly flagged.
- [ ] Stage 6 never modifies the Stage-5 trajectory.
- [ ] Stage 6 never fits the TLE to the observations.

## Engineering

- [ ] Stage-1 through Stage-6 tests pass.
- [ ] Registry tests pass.
- [ ] No network required for unit tests.
- [ ] Public functions have docstrings.
- [ ] Logging is used consistently.

---

# 46. Do Not Implement Yet

Do not implement in Stage 6:

```text
automatic anomaly classification
station-specific bias attribution
clock-error inference
calibration correction
RFI attribution
satellite manoeuvre claims
sensor-health scoring
multi-satellite population statistics
```

These belong to Stage 7 or later analyses.

---

# 47. Expected End State

The public pipeline becomes:

```bash
tart-stations-update ...
tart-stations-list ...
tart-common-gnss ...
tart-retrieve ...
tart-localise ...
tart-triangulate ...
tart-track ...
tart-validate ...
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
3-D measured positions + covariance
      |
tart-track
      |
TLE-independent smoothed trajectory
      |
tart-validate
      |
observed - TLE residuals
Cartesian + RTN residual statistics
measurement uncertainty
reference-quality diagnostics
      |
Stage 7
      |
systematic / station / anomaly interpretation
```

---

# 48. Final Deliverables

The coding agent must deliver:

```text
src/tart_gnss/stage6/*
tests/stage6/*
docs/stage6_validation.md
README.md update
pyproject.toml entry-point/dependency updates
tests/test_cli_names.py update
tart_gnss_stage6_implementation.md
```

Also provide:

```text
synthetic Stage-6 validation example
real-data Stage-6 example when an appropriate TLE is available
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
known injected RTN offset
recovered radial residual
recovered along-track residual
recovered cross-track residual
known velocity offset
recovered RTN velocity residual
3-D residual norm
angular residual
measurement-conditional Mahalanobis
reference-aware Mahalanobis if ensemble enabled
TLE epoch
TLE age
precision-validation status
```

The implementation is complete only when the selected TLE is introduced **after** the Stage-5 TLE-independent trajectory has been produced, and the resulting residuals remain scientifically interpretable without treating the TLE as exact ground truth.
