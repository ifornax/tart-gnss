# TART-GNSS Stage 3 Implementation Specification
## Visibility-Domain Direction Estimation and Image-Domain Validation

**Repository:** `git@github.com:ifornax/tart-gnss.git`  
**Target branch:** create a new implementation branch from `main`  
**Stage status before this work:** Stage 1 and Stage 2 are complete  
**This specification implements:** Stage 3 only

---

# 1. Objective

Implement **Stage 3: satellite direction estimation from real TART interferometric visibility data**.

Stage 3 must estimate the apparent direction of a GNSS satellite independently at each TART station using the complex interferometric visibilities produced/retrieved by Stage 2.

The primary estimator must operate directly in the **visibility domain**.

An image-domain estimator must also be implemented as an independent validation path.

Stage 3 must produce:

1. time-resolved satellite direction estimates;
2. angular uncertainty estimates;
3. goodness-of-fit and data-quality metrics;
4. image-domain cross-checks;
5. reproducible machine-readable outputs for Stage 4;
6. diagnostic plots suitable for scientific validation.

Do **not** combine geographically separated TART stations into a single interferometer.

Each TART station must be processed independently.

Stage 3 must consume the existing Stage-2 `manifest.json` and the Measurement Sets/native files already produced by the pipeline.

---

# 2. Scientific Principle

For one dominant point-like satellite, the complex visibility measured by antennas \(p\) and \(q\) can be approximated as

\[
V_{pq}(t,\nu)
=
A(t,\nu)
\exp
\left[
-\frac{2\pi i}{\lambda}
\mathbf{b}_{pq}(t)\cdot\mathbf{s}(t)
\right]
+
\epsilon_{pq}(t,\nu),
\]

where

- \(V_{pq}\) is the observed complex visibility;
- \(A\) is the apparent satellite signal amplitude;
- \(\lambda=c/\nu\);
- \(\mathbf{b}_{pq}\) is the baseline vector;
- \(\mathbf{s}\) is the source direction unit vector;
- \(\epsilon_{pq}\) represents thermal noise, calibration residuals, RFI, and modelling error.

For a local tangent-plane representation,

\[
\mathbf{s}(l,m)
=
\begin{bmatrix}
l \\
m \\
\sqrt{1-l^2-m^2}
\end{bmatrix},
\]

subject to

\[
l^2+m^2 \le 1.
\]

The visibility-domain estimator must infer

\[
\hat l,\hat m
\]

and transform those estimates to the station-local

\[
\hat{\mathrm{Az}}, \hat{\mathrm{El}}.
\]

The estimated direction must not be forced to the TLE direction.

TLE predictions may be used for:

- selecting the target satellite;
- defining a broad search region;
- source association;
- plotting comparisons after the measurement estimate has been obtained.

They must **not** be used as a strong prior that forces the fitted position toward the predicted direction.

---

# 3. Stage Boundaries

Existing pipeline:

```text
Stage 1
    Common GNSS visibility prediction
          |
          v
Stage 2
    Retrieve + validate + synchronise observations
    + one Measurement Set per TART station
          |
          v
Stage 3   <-- IMPLEMENT HERE
    Direction estimation from visibilities
          |
          v
Stage 4
    Multi-station 3-D localisation
```

Stage 3 must not implement:

- 3-D triangulation;
- multi-station position inference;
- trajectory filtering;
- orbital fitting;
- TLE residual science analysis.

Those belong to later stages.

---

# 4. Repository Integration

Use the current package structure.

Add a new module tree:

```text
src/tart_gnss/stage3/
├── __init__.py
├── cli.py
├── loader.py
├── geometry.py
├── visibility_model.py
├── estimator.py
├── imaging.py
├── uncertainty.py
├── quality.py
├── output.py
├── plotting.py
└── errors.py
```

Add tests under:

```text
tests/stage3/
├── test_loader.py
├── test_geometry.py
├── test_visibility_model.py
├── test_estimator.py
├── test_uncertainty.py
├── test_quality.py
├── test_output.py
├── test_cli.py
└── test_synthetic_recovery.py
```

Do not rewrite Stage-1 or Stage-2 modules unless a minimal compatibility change is absolutely necessary.

Any required compatibility change must be:

- backward compatible;
- documented;
- covered by tests.

---

# 5. CLI

Implement a new executable:

```bash
tart-localise
```

Suggested usage:

```bash
tart-localise \
    --manifest data/observations/20260625/NORAD_40534/manifest.json \
    --station mu-udm \
    --output-dir results/stage3 \
    --method visibility
```

Support:

```bash
tart-localise \
    --manifest <manifest.json> \
    --all-stations \
    --output-dir results/stage3
```

Required options:

```text
--manifest
--station
--all-stations
--output-dir
--method
--time-bin
--min-baselines
--search-radius-deg
--image-validation
--verbose
```

Allowed methods initially:

```text
visibility
image
both
```

Default:

```text
--method visibility
```

`--image-validation` should run the image-domain estimator even if the primary result is visibility-domain.

---

# 6. Stage-2 Manifest Input

Stage 3 must treat Stage-2 `manifest.json` as the authoritative input.

Validate before processing:

- `three_station_valid`;
- station identity;
- Measurement Set path;
- observation frequency;
- time coverage;
- calibration state;
- number of antennas;
- flagged/zero-gain antennas;
- predicted target satellite;
- synchronization metadata.

If the full three-station Stage-2 run is invalid, Stage 3 may still allow processing of an individually valid station only when explicitly requested.

In that case, clearly mark:

```json
"stage2_three_station_valid": false
```

in Stage-3 outputs.

Never silently reinterpret an invalid Stage-2 dataset as valid.

---

# 7. Measurement Set Loading

Implement a robust MS reader.

Read at minimum:

```text
TIME
DATA
UVW
ANTENNA1
ANTENNA2
FLAG
FLAG_ROW
WEIGHT
SIGMA
```

and relevant subtables:

```text
ANTENNA
SPECTRAL_WINDOW
FIELD
OBSERVATION
```

Use `python-casacore` if already available through the existing Stage-2 dependency path.

The loader must return a clean internal representation, preferably using NumPy arrays or a lightweight dataclass.

Example:

```python
@dataclass
class VisibilityDataset:
    time_mjd_seconds: np.ndarray
    antenna1: np.ndarray
    antenna2: np.ndarray
    uvw_m: np.ndarray
    vis: np.ndarray
    flag: np.ndarray
    weight: np.ndarray
    frequency_hz: np.ndarray
    station_id: str
```

Avoid introducing a heavy dependency unless necessary.

---

# 8. Time Binning

Satellite positions vary continuously, so direction estimation must be time-resolved.

Implement configurable temporal binning.

Default:

```text
--time-bin 1.0
```

seconds, but do not assume the native cadence is exactly 1 s.

For each bin:

1. select visibility rows falling inside the interval;
2. remove flagged rows;
3. reject invalid baselines;
4. estimate one apparent source direction;
5. compute uncertainty and quality metrics.

The output time must be the weighted or arithmetic midpoint of the measurements actually used.

Do not fabricate timestamps on a fixed grid if no measurements exist.

---

# 9. Baseline Selection

Exclude:

- auto-correlations;
- fully flagged baselines;
- baselines involving invalid/zero-gain antennas when identified by Stage 2;
- rows with NaN or infinite values;
- rows with zero or invalid weights.

Record:

```text
n_baselines_total
n_baselines_used
fraction_flagged
```

for each time bin.

Add configurable:

```text
--min-baselines
```

If fewer than the required number remain, no scientific direction estimate should be reported.

Instead output a quality flag.

---

# 10. Coordinate Geometry

Implement conversions between:

```text
(l, m)
unit vector
azimuth/elevation
ENU
ITRS/ECEF where required
```

Be explicit about coordinate conventions.

Recommended local frame:

```text
x = East
y = North
z = Up
```

Azimuth convention:

```text
0° = North
90° = East
```

Elevation:

```text
0° = horizon
90° = zenith
```

Add unit tests for known cardinal directions.

All coordinate transformations must be documented mathematically.

---

# 11. Visibility Forward Model

Implement:

```python
model_visibility(...)
```

that predicts complex visibilities for an assumed direction.

For narrow-band first implementation:

\[
V_k^{\rm model}
=
A\,
e^{-2\pi i \mathbf{u}_k\cdot\mathbf{s}},
\]

where baseline coordinates are in wavelengths.

Amplitude may initially be modelled as a nuisance complex scalar:

\[
C = a + ib.
\]

For a given candidate direction, solve analytically or numerically for the best complex amplitude.

Recommended formulation:

\[
V_k \approx C\,e^{-i\phi_k(l,m)}.
\]

For fixed \(l,m\), estimate

\[
\hat C
=
\frac{\sum_k w_k V_k e^{+i\phi_k}}
{\sum_k w_k}.
\]

Then minimise the weighted residual:

\[
\chi^2(l,m)
=
\sum_k
w_k
\left|
V_k
-
\hat C e^{-i\phi_k(l,m)}
\right|^2.
\]

This eliminates two nuisance parameters and leaves a two-dimensional optimisation problem.

---

# 12. Primary Visibility-Domain Estimator

Implement a robust two-stage optimisation.

## 12.1 Coarse search

Search a bounded \(l,m\) region.

Use the TLE-predicted direction only to define an optional broad search centre.

Default search radius:

```text
--search-radius-deg 10
```

The result must still be determined by the visibilities.

If a reliable blind full-sky search is practical, support it as:

```text
--blind
```

but this is not required for the first implementation.

## 12.2 Local refinement

Refine the best coarse-grid candidate with:

```text
scipy.optimize.minimize
```

or `least_squares`.

Recommended methods:

```text
L-BFGS-B
Powell
least_squares
```

Use bounded parameters.

Do not allow:

\[
l^2+m^2 > 1.
\]

Return:

```text
l_hat
m_hat
az_hat_deg
el_hat_deg
complex_amplitude_hat
chi2
reduced_chi2
success
optimizer_status
```

---

# 13. Robust Estimation

Real TART data may contain bad baselines or unmodelled sources.

Support robust loss.

Minimum required options:

```text
linear
huber
soft_l1
```

Default:

```text
huber
```

If using a residual formulation with real and imaginary components:

\[
r =
[
\Re(V-V_m),
\Im(V-V_m)
].
\]

Use robust loss through `scipy.optimize.least_squares`.

Record the chosen loss function in provenance.

---

# 14. Phase-Wrapping

Do not fit baseline phase differences by naïvely subtracting wrapped phases.

Fit the complex visibility directly.

This avoids phase-wrap discontinuities.

Any phase-only diagnostic may use wrapped residuals for plotting, but not as the primary objective function.

---

# 15. Multi-Frequency Handling

TART data are narrow-band, but the implementation should not hard-code a single frequency channel.

If multiple channels exist:

\[
\phi_{k,\nu}
=
2\pi
\frac{\nu}{c}
\mathbf{b}_k\cdot\mathbf{s}.
\]

The model must use each channel's actual frequency.

Allow averaging only if:

- phase smearing is negligible;
- averaging is explicitly requested;
- the operation is recorded.

---

# 16. Near-Field Check

GNSS satellites are not mathematically at infinite distance.

Before finalising the far-field formulation, implement a numerical diagnostic comparing the array's maximum baseline \(B\) and satellite range \(R\).

Compute a Fresnel-style far-field metric such as

\[
\eta = \frac{B^2}{\lambda R}.
\]

Record this diagnostic.

Do not silently assume the plane-wave approximation without checking it.

If the expected near-field phase correction is negligible compared with the measurement precision, document that result.

If not negligible, add a spherical-wave model option:

\[
V_{pq}
\propto
\exp
\left[
-\frac{2\pi i}{\lambda}
\left(
|\mathbf{x}-\mathbf{r}_p|
-
|\mathbf{x}-\mathbf{r}_q|
\right)
\right].
\]

For Stage 3, only implement this second model if the diagnostic demonstrates that it is needed.

Do not over-engineer it without evidence.

---

# 17. Direction Uncertainty

Every successful fit must return uncertainty.

Implement at least one defensible frequentist covariance estimate.

Recommended first method:

\[
\Sigma_{lm}
\approx
(J^\mathsf{T} W J)^{-1}\hat\sigma^2.
\]

Where:

- \(J\) is the Jacobian of the residuals;
- \(W\) contains visibility weights;
- \(\hat\sigma^2\) comes from residual variance.

Return:

```text
sigma_l
sigma_m
cov_lm
```

Then propagate uncertainty into azimuth/elevation numerically or analytically.

Return:

```text
sigma_az_deg
sigma_el_deg
cov_az_el
```

Be cautious near zenith, where azimuth becomes poorly conditioned.

Flag such cases rather than reporting misleadingly small azimuth uncertainty.

---

# 18. Bootstrap Validation

Implement an optional baseline bootstrap:

```text
--bootstrap N
```

Default:

```text
0
```

When enabled:

1. resample baselines with replacement;
2. repeat localisation;
3. estimate empirical confidence intervals.

Output at minimum:

```text
az_p16
az_p50
az_p84
el_p16
el_p50
el_p84
```

This provides a useful check on the local-Hessian covariance.

---

# 19. Image-Domain Validation

Implement an independent image-based direction estimator.

This is a validation estimator, not the default scientific estimator.

Possible route:

1. form dirty image from calibrated visibilities;
2. search for the dominant source near the expected GNSS region;
3. fit the image peak;
4. convert image coordinates to azimuth/elevation.

Do not require CLEAN for the first implementation unless demonstrably necessary.

Produce:

```text
az_image_deg
el_image_deg
peak_value
image_snr
```

Then calculate estimator disagreement:

\[
\Delta\theta_{\rm vis-img}
=
\cos^{-1}
\left(
\mathbf{s}_{\rm vis}\cdot\mathbf{s}_{\rm img}
\right).
\]

This is a Stage-3 internal consistency check.

Do not compare against the TLE as the only estimator validation.

---

# 20. Signal-to-Noise Metric

Define and document an operational SNR.

For example:

\[
{\rm SNR}
=
\frac{|\hat C|}
{\sigma_{\rm residual}}.
\]

Or derive an equivalent robust residual statistic.

Do not report an undefined or ad-hoc quantity labelled "SNR".

Store the exact definition in code documentation.

---

# 21. Goodness of Fit

For every estimate compute at minimum:

```text
chi2
reduced_chi2
rms_complex_residual
median_abs_residual
max_abs_residual
```

If the statistical assumptions required for formal \(\chi^2\) interpretation are not satisfied, still use the statistic as a diagnostic but do not overstate its probabilistic meaning.

---

# 22. Quality Flags

Implement deterministic quality flags.

Suggested values:

```text
GOOD
LOW_BASELINE_COUNT
HIGH_FLAGGED_FRACTION
FIT_FAILED
BOUNDARY_SOLUTION
LOW_SNR
HIGH_RESIDUAL
IMAGE_DISAGREEMENT
ZENITH_AZ_UNCERTAIN
INVALID_GEOMETRY
CALIBRATION_WARNING
```

Multiple flags should be supported.

Do not collapse all problems into a single Boolean.

Output:

```json
"quality_flags": []
```

---

# 23. Stage-3 Output Schema

Create one result directory per target run.

Recommended:

```text
results/stage3/
└── 20260625/
    └── NORAD_40534/
        ├── stage3_manifest.json
        ├── directions.parquet
        ├── directions.csv
        ├── diagnostics/
        │   ├── mu-udm_track.png
        │   ├── mu-udm_fit_quality.png
        │   ├── mu-udm_vis_image_offset.png
        │   └── ...
        └── logs/
            └── stage3.log
```

Primary tabular output fields:

```text
time_utc
time_mjd
station_id
norad_id

l
m
sigma_l
sigma_m
cov_lm

az_obs_deg
el_obs_deg
sigma_az_deg
sigma_el_deg
cov_az_el

complex_amp_real
complex_amp_imag
apparent_amplitude

chi2
reduced_chi2
rms_complex_residual
median_abs_residual
snr

n_baselines_total
n_baselines_used
flagged_fraction

az_image_deg
el_image_deg
image_snr
vis_image_separation_deg

fit_success
optimizer_status
quality_flags
```

TLE-related columns may be appended for convenience:

```text
az_tle_deg
el_tle_deg
```

but these must be marked as **reference predictions**, not fitted quantities.

Do not yet make TLE residuals the primary Stage-3 science output.

---

# 24. Stage-3 Manifest

Write:

```text
stage3_manifest.json
```

containing:

```json
{
  "schema_version": 1,
  "stage": 3,
  "input_manifest": "...",
  "input_manifest_sha256": "...",
  "norad_id": "...",
  "stations": {},
  "estimator": {},
  "time_binning": {},
  "quality_thresholds": {},
  "software": {},
  "provenance": {},
  "outputs": {}
}
```

Include:

- git commit;
- package version;
- Python version;
- NumPy/SciPy versions;
- casacore version;
- command line;
- timestamp;
- input hashes;
- estimator options;
- robust loss;
- search range;
- time bin;
- uncertainty method;
- bootstrap count.

Follow the same provenance philosophy already used in Stage 1 and Stage 2.

---

# 25. Diagnostic Plots

Implement at minimum:

## 25.1 Measured direction vs time

Two panels:

```text
Azimuth vs time
Elevation vs time
```

Plot:

- visibility-domain estimate;
- optional image-domain estimate;
- TLE prediction as a reference line;
- uncertainty interval where visible.

## 25.2 Sky track

Plot:

```text
azimuth vs elevation
```

or a suitable local-sky projection.

Show:

- visibility estimate;
- image estimate;
- TLE prediction.

## 25.3 Visibility-image disagreement

Plot:

\[
\Delta\theta_{\rm vis-img}(t).
\]

## 25.4 Fit quality

Plot:

```text
SNR
reduced chi2
baseline count
flagged fraction
```

against time.

Use matplotlib only.

Plots must be generated from saved Stage-3 tables, not only from transient in-memory objects.

---

# 26. Calibration Metadata

Carry Stage-2 calibration information forward.

For every station, record:

- calibration applied/unapplied;
- zero-gain antennas;
- calibration changes;
- Stage-2 validation warnings.

If the station has significant calibration issues, add:

```text
CALIBRATION_WARNING
```

to relevant Stage-3 results.

Do not attempt a new calibration pipeline in Stage 3.

That is a separate scientific task.

---

# 27. Timing Metadata

Stage 2 already establishes that separate TART stations are not perfectly simultaneous.

Stage 3 estimates directions using the timestamps belonging to each station independently.

Do not force station measurements onto identical timestamps in Stage 3.

Store the actual estimate time.

Time alignment belongs to Stage 4.

---

# 28. TLE Use

Stage 3 must maintain the distinction:

```text
prediction != measurement
```

TLEs may:

- identify the expected GNSS satellite;
- define a broad search region;
- assist in plotting;
- reject obviously unrelated detections after localisation.

TLEs must not:

- set the fitted result;
- act as an artificially tight prior;
- replace a failed visibility fit;
- be copied into output as an "observed" position.

---

# 29. Synthetic Validation

Create a synthetic visibility generator for tests.

Given:

```text
array geometry
frequency
source azimuth/elevation
complex amplitude
noise level
flags
```

generate visibilities using the same forward model.

Required tests:

### Noise-free recovery

Recover input direction to numerical precision.

### Gaussian-noise recovery

Recovered direction should be statistically consistent with the injected value.

### Flagging

Direction should remain recoverable after randomly flagging a reasonable subset of baselines.

### Bad baselines

Robust loss should outperform ordinary least squares when outliers are injected.

### Boundary case

Source near horizon should not produce invalid \(l,m\).

### Zenith

Azimuth uncertainty must be handled appropriately.

### TLE-offset test

Inject a source deliberately offset from the search-centre/TLE position and verify that the estimator recovers the injected position rather than returning the TLE position.

This last test is mandatory.

---

# 30. Real-Data Acceptance Dataset

Use an existing valid Stage-2 dataset, preferably the existing:

```text
2026-06-25
NORAD 40534
```

if the local files are available.

The existing Stage-2 report indicates:

- three valid station Measurement Sets;
- 24 antennas per station before flagging;
- real common overlap approximately 645 s;
- Rhodes has significant zero-gain antenna flagging;
- Mauritius and Madagascar show gain-solution changes.

This makes it a useful first real-data stress test.

Do not commit large Measurement Sets to git.

Keep real-data tests optional/skipped when files are absent.

---

# 31. Acceptance Criteria

Stage 3 is complete only when all of the following are true.

## Functional

- [ ] `tart-localise` exists.
- [ ] It reads a Stage-2 manifest.
- [ ] It processes one station.
- [ ] It processes all available stations independently.
- [ ] It estimates direction directly from complex visibilities.
- [ ] It returns uncertainties.
- [ ] It produces image-domain cross-checks.
- [ ] It writes CSV/Parquet output.
- [ ] It writes `stage3_manifest.json`.
- [ ] It generates diagnostic figures.

## Scientific

- [ ] Synthetic sources are recovered accurately.
- [ ] TLE search-centre offsets do not bias the fitted solution.
- [ ] Complex visibilities are fitted directly.
- [ ] Phase wrapping is not used as the primary residual.
- [ ] Flagged/invalid baselines are excluded.
- [ ] Visibility estimator and image estimator broadly agree for clean synthetic tests.
- [ ] Uncertainty estimates are numerically stable.
- [ ] Near-field validity is explicitly checked and recorded.

## Engineering

- [ ] Existing Stage-1 tests pass.
- [ ] Existing Stage-2 tests pass.
- [ ] Registry tests pass.
- [ ] New Stage-3 tests pass.
- [ ] No network access is required for unit tests.
- [ ] Large observational data are not added to git.
- [ ] Code is typed where practical.
- [ ] Public functions have docstrings.
- [ ] Logging uses the existing project logging approach.

---

# 32. API Design

Prefer small, testable functions.

Example interfaces:

```python
def load_station_ms(
    ms_path: Path,
    station_id: str,
) -> VisibilityDataset:
    """Load and validate one station Measurement Set."""
```

```python
def predict_visibilities(
    uvw_m: np.ndarray,
    frequency_hz: np.ndarray,
    l: float,
    m: float,
    complex_amplitude: complex,
) -> np.ndarray:
    """Predict narrow-field complex visibilities for a point source."""
```

```python
def estimate_direction(
    dataset: VisibilityDataset,
    initial_az_deg: float | None = None,
    initial_el_deg: float | None = None,
    search_radius_deg: float = 10.0,
    robust_loss: str = "huber",
) -> DirectionEstimate:
    """Estimate source direction from complex interferometric visibilities."""
```

```python
def estimate_direction_uncertainty(
    result: DirectionEstimate,
    dataset: VisibilityDataset,
) -> DirectionCovariance:
    """Estimate local covariance of the fitted source direction."""
```

```python
def image_direction(
    dataset: VisibilityDataset,
    search_region: SearchRegion | None = None,
) -> ImageDirectionEstimate:
    """Estimate source direction from a dirty-image peak."""
```

---

# 33. Dataclasses

Suggested:

```python
@dataclass
class DirectionEstimate:
    time_utc: datetime
    station_id: str
    norad_id: int
    l: float
    m: float
    az_deg: float
    el_deg: float
    complex_amplitude: complex
    chi2: float
    reduced_chi2: float
    snr: float
    n_baselines_used: int
    success: bool
    quality_flags: list[str]
```

and

```python
@dataclass
class DirectionCovariance:
    sigma_l: float
    sigma_m: float
    cov_lm: float
    sigma_az_deg: float
    sigma_el_deg: float
    cov_az_el: float
```

Keep the output schema stable and versioned.

---

# 34. Dependencies

Prefer existing dependencies plus:

```text
scipy
pandas
pyarrow
matplotlib
```

only if not already present.

Do not add:

- TensorFlow;
- PyTorch;
- JAX;
- probabilistic-programming frameworks;

for Stage 3.

The primary fit does not require them.

---

# 35. Numerical Considerations

Be careful about:

- radians vs degrees;
- metres vs wavelengths;
- channel frequencies;
- UVW convention;
- complex conjugation convention;
- baseline orientation;
- source direction sign convention.

Create explicit tests capable of detecting a sign flip.

A source injected East of phase centre must not be recovered West because of an unnoticed Fourier convention mismatch.

---

# 36. UVW Convention Validation

Do not assume the Measurement Set UVW convention without testing it.

Create a deterministic synthetic/reference test that confirms:

```text
visibility phase sign
UVW orientation
(l,m) sign
```

against the adopted measurement equation.

Document the convention in:

```text
docs/stage3_direction_estimation.md
```

This is mandatory.

---

# 37. Documentation

Add:

```text
docs/stage3_direction_estimation.md
```

It must include:

1. scientific objective;
2. visibility measurement equation;
3. coordinate conventions;
4. estimator derivation;
5. complex-amplitude elimination;
6. robust loss;
7. uncertainty calculation;
8. imaging validation;
9. quality flags;
10. near-field diagnostic;
11. output schema;
12. command examples;
13. limitations.

Update the main `README.md` roadmap:

```text
Stage 3 — implemented
```

only after acceptance criteria are met.

---

# 38. Reproducibility

Every Stage-3 run must be reproducible from:

```text
Stage-2 manifest
Measurement Sets
Stage-3 manifest
git commit
command line
configuration
```

Do not rely on hidden runtime state.

If randomness is used for bootstrap analysis, record:

```text
random_seed
```

in the manifest.

---

# 39. Logging

Use `logging`.

At INFO level report:

```text
station being processed
input MS
number of rows
time range
frequency range
number of valid antennas
number of valid baselines
number of time bins
fit success fraction
median SNR
median visibility-image offset
output paths
```

At DEBUG level allow detailed optimizer diagnostics.

Do not print large visibility arrays.

---

# 40. Failure Behaviour

Use explicit exceptions in:

```text
stage3/errors.py
```

Suggested:

```text
Stage3Error
InvalidManifestError
MissingMeasurementSetError
InvalidVisibilityDataError
InsufficientBaselinesError
DirectionFitError
CoordinateTransformError
```

CLI exit codes should distinguish:

```text
0  successful scientific run
2  invalid input/configuration
3  completed but no valid scientific estimates
1  unexpected execution failure
```

Match Stage-2 philosophy where practical.

---

# 41. Do Not Implement Yet

Do not implement in this task:

- multi-station triangulation;
- satellite ECEF/ITRS position recovery;
- Bayesian \(x,y,z\) inference;
- Kalman filtering;
- orbital dynamics;
- along-track/cross-track residuals;
- anomaly classification;
- station-bias inference;
- automatic calibration correction.

These belong to Stages 4–7.

---

# 42. Expected End State

After implementation, the following should work:

```bash
tart-localise \
  --manifest data/observations/20260625/NORAD_40534/manifest.json \
  --all-stations \
  --method visibility \
  --image-validation \
  --time-bin 1 \
  --output-dir results/stage3
```

and produce three independent station tracks:

```text
Mauritius:
    Az_obs(t), El_obs(t), covariance

Madagascar:
    Az_obs(t), El_obs(t), covariance

South Africa:
    Az_obs(t), El_obs(t), covariance
```

These become the direct input to Stage 4:

```text
station positions
+
measured line-of-sight directions
+
direction covariance
        |
        v
3-D GNSS satellite localisation
```

---

# 43. Final Deliverables

The coding agent must deliver:

```text
src/tart_gnss/stage3/*
tests/stage3/*
docs/stage3_direction_estimation.md
README.md update
pyproject.toml dependency/entry-point updates
```

plus one small example output generated from synthetic data.

Do not commit large real TART datasets.

Commit messages should separate logically distinct work where practical.

Before completion run:

```bash
pytest
```

and report:

```text
number passed
number failed
number skipped
```

Also run a synthetic end-to-end demonstration and report:

```text
injected az/el
recovered az/el
angular error
estimated uncertainty
image-domain result
visibility-image disagreement
```

The implementation is complete only when the recovered synthetic direction is demonstrably measurement-driven and does not collapse to the TLE/search-centre position.
