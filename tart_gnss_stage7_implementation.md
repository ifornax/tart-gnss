# TART-GNSS Stage 7 Implementation Specification
## Sensor-Network Diagnosis, Systematic-Effect Attribution, and Anomaly Evidence

**Repository:** `git@github.com:ifornax/tart-gnss.git`  
**Base branch:** `main`  
**Stage status before this work:** Stages 1–6 implemented  
**This specification implements:** Stage 7 only  
**Canonical Stage-7 CLI:** `tart-diagnose`

---

# 0. Permanent CLI Naming Convention

Every public CLI in this repository must begin with:

```text
tart-
```

Current canonical commands:

```text
tart-stations-update
tart-stations-list
tart-common-gnss
tart-retrieve
tart-localise
tart-triangulate
tart-track
tart-validate
```

Stage 7 must add:

```text
tart-diagnose
```

Do not add public commands with any other naming pattern.

Update `pyproject.toml`, `README.md`, docs, tests, and examples accordingly. Historical provenance files must remain unchanged.

---

# 1. Objective

Implement **Stage 7: sensor-network diagnosis and systematic-effect attribution**.

Stages 3–6 answer where the satellite is, how uncertain that estimate is, and how the measured trajectory differs from the TLE prediction.

Stage 7 must answer:

```text
What evidence suggests where the discrepancy comes from?
```

The stage must combine information from:

```text
Stage 2  data/calibration/timing metadata
Stage 3  station-level direction fits and quality
Stage 4  triangulation geometry and station line residuals
Stage 5  trajectory innovations and dynamics diagnostics
Stage 6  observed-minus-TLE Cartesian and RTN residuals
```

The goal is to distinguish signatures consistent with station-specific bias, calibration problems, timing mismatch, poor network geometry, localisation failure, trajectory-model mismatch, TLE/reference uncertainty, or a network-consistent unexplained residual.

Stage 7 must not automatically claim physical cause.

---

# 2. Scientific Principle

A final residual may be considered as a mixture of:

```text
station effects
measurement effects
geometry effects
trajectory-model effects
TLE/reference effects
unexplained remainder
```

Stage 7 must therefore examine each layer separately.

The implementation should favour:

```text
diagnostic evidence
```

over:

```text
automatic causal attribution
```

---

# 3. Evidence Hierarchy

Implement the following evidence levels:

```text
Level 0 — NORMAL
No statistically meaningful anomaly after quality checks.

Level 1 — STATISTICAL_ANOMALY
A residual or innovation is unusual relative to its stated uncertainty.

Level 2 — STATION_CORRELATED
The anomaly is associated with one station or disappears when that station is removed.

Level 3 — INSTRUMENT_CORRELATED
The anomaly correlates with calibration, flagged antennas, gain changes, timing offsets,
image/visibility disagreement, or another instrument/data-quality diagnostic.

Level 4 — NETWORK_CONSISTENT
The anomaly persists across station subsets and is not explained by known station/instrument effects.

Level 5 — UNEXPLAINED
The anomaly is network-consistent, statistically strong, persistent, and remains unexplained
by available diagnostics and TLE/reference-quality indicators.
```

These are evidence levels, not physical-cause labels.

Do not equate Level 5 with a manoeuvre or any other physical event.

---

# 4. Stage Boundary

```text
Stage 1  common visibility
Stage 2  data retrieval and validation
Stage 3  direction estimation
Stage 4  3-D localisation
Stage 5  trajectory inference
Stage 6  TLE validation
Stage 7  systematic/network diagnosis
```

Stage 7 must not modify upstream scientific outputs.

---

# 5. Required Inputs

Primary input:

```text
stage6_manifest.json
```

Stage 7 must trace provenance backward and load associated outputs from:

```text
Stage 6: residuals.parquet, summary.json
Stage 5: trajectory.parquet, innovations.parquet, stage5_manifest.json
Stage 4: positions.parquet, station_contributions.parquet, stage4_manifest.json
Stage 3: directions.parquet, stage3_manifest.json
Stage 2: manifest.json
```

Use recorded provenance links and hashes wherever available.

Do not discover unrelated files heuristically when a manifest path exists.

---

# 6. CLI

Implement:

```bash
tart-diagnose
```

Example:

```bash
tart-diagnose \
    --stage6-manifest results/stage6/20260625/NORAD_40534/stage6_manifest.json \
    --output-dir results/stage7
```

Required/supported options:

```text
--stage6-manifest
--output-dir
--quality-policy
--leave-one-out
--pairwise
--correlation-method
--min-samples
--anomaly-threshold
--min-persistent-epochs
--no-plots
--verbose
```

Defaults:

```text
--quality-policy strict
--leave-one-out true
--pairwise true
--correlation-method spearman
```

---

# 7. Repository Integration

Add:

```text
src/tart_gnss/stage7/
├── __init__.py
├── cli.py
├── loader.py
├── alignment.py
├── influence.py
├── correlations.py
├── attribution.py
├── scoring.py
├── quality.py
├── statistics.py
├── output.py
├── plotting.py
└── errors.py
```

Tests:

```text
tests/stage7/
├── __init__.py
├── test_loader.py
├── test_alignment.py
├── test_influence.py
├── test_correlations.py
├── test_attribution.py
├── test_scoring.py
├── test_quality.py
├── test_statistics.py
├── test_output.py
├── test_cli.py
└── test_synthetic_diagnosis.py
```

Documentation:

```text
docs/stage7_diagnosis.md
```

---

# 8. Time Alignment Across Stages

Different stage outputs may have different row counts and timestamps.

Do not join by row number.

Use:

```text
time_utc
station_id
norad_id
```

with configurable tolerance where necessary.

For each Stage-6 epoch, link the Stage-6 residual, Stage-5 innovation/trajectory diagnostics, Stage-4 triangulation diagnostics, Stage-3 station measurements, and Stage-2 calibration/timing metadata.

Record missing provenance explicitly.

---

# 9. Leave-One-Station-Out Analysis

This is a primary Stage-7 diagnostic.

For each station i, compute a leave-one-station-out solution:

```text
x_hat_minus_i(t)
```

Prefer reusing Stage-4 triangulation functions directly.

Do not duplicate triangulation mathematics.

Compute the difference between the full solution and the leave-one-out solution.

Output:

```text
loo_shift_3d_m
loo_shift_radial_m
loo_shift_along_track_m
loo_shift_cross_track_m
```

when a Stage-6 RTN basis is available.

A station with consistently large leave-one-out influence is diagnostically important.

---

# 10. Three-Station Network Caveat

With exactly three stations, leaving one out leaves two stations.

Two-station localisation is geometrically weaker.

Therefore leave-one-out results must include:

```text
n_stations_remaining
geometry_condition_number
crossing_angle
degraded_geometry
```

Do not interpret a large leave-one-out shift without accounting for geometry degradation.

---

# 11. Normalised Station Influence

Define a normalised influence statistic using the leave-one-out covariance where available.

Output:

```text
loo_mahalanobis
```

Do not call it a probability unless assumptions are valid.

If covariance is unavailable or singular, report only the physical shift.

---

# 12. Pairwise Station Diagnostics

For networks with at least three stations, compute pairwise two-station diagnostic localisations.

For each station pair, estimate the closest-approach position or equivalent two-line solution.

Output:

```text
pair_station_1
pair_station_2
pair_position
pair_line_separation_m
pair_geometry_angle_deg
```

These are diagnostic-only solutions.

---

# 13. Station Direction Residual Diagnostics

For each station, calculate its measured direction relative to the final Stage-5/Stage-6 trajectory projected back to that station.

Compute:

```text
direction_residual_az_deg
direction_residual_el_deg
direction_residual_angle_deg
```

using the Stage-5 smoothed trajectory as the network-derived reference.

This helps isolate station-specific direction bias.

---

# 14. Station Direction Bias

Estimate robust station-level bias:

```text
median az residual
median el residual
median angular residual
MAD scale
```

Optionally fit simple descriptive models versus:

```text
elevation
azimuth
time
```

Do not automatically calibrate these biases away.

---

# 15. Calibration Correlation Analysis

From Stage 2 / Stage 3 gather:

```text
zero-gain antenna count
flagged fraction
gain-solution changes
calibration warnings
Stage-3 SNR
Stage-3 reduced chi-square
Stage-3 image/visibility disagreement
baseline count
```

Test correlations against:

```text
station direction residual
Stage-4 line residual
leave-one-out influence
Stage-6 RTN residual
```

Default:

```text
Spearman rank correlation
```

Output:

```text
correlation coefficient
p-value
n_samples
```

Treat p-values as descriptive evidence, not automated truth.

---

# 16. Timing Correlation Analysis

Use:

```text
Stage-4 station_time_offset_s
Stage-4 TIME_MISMATCH flags
Stage-2 cross-station timing diagnostics
```

Compare with:

```text
along-track residual
station line residual
direction residual
```

Timing-like effects should preferentially appear in along-track behaviour.

Do not automatically label an along-track residual as a clock error.

---

# 17. Geometry Dependence

Analyse whether residuals increase with poor triangulation geometry.

Inputs:

```text
condition_number
minimum crossing angle
maximum crossing angle
number of stations
```

Compare with:

```text
Stage-4 position uncertainty
Stage-6 residual norm
RTN residuals
leave-one-out shift
```

Expected pattern:

```text
poorer geometry -> larger uncertainty
```

A residual that scales with geometry is weaker evidence for a physical anomaly.

---

# 18. Stage-5 Dynamics Diagnostics

Use:

```text
NIS
innovation outliers
DYNAMICS_MISMATCH
HIGH_NIS
FILTER_DIVERGENCE
innovation autocorrelation
```

Compare with Stage-6 residuals.

If a Stage-6 anomaly occurs at the same time as Stage-5 dynamics mismatch, classify it as:

```text
TRAJECTORY_MODEL_CORRELATED
```

rather than immediately network-consistent.

---

# 19. TLE Reference Diagnostics

Use Stage-6 fields:

```text
tle_age_days
precision_validation
TLE_STALE
TLE_BACKPROPAGATED
reference ensemble dispersion
REFERENCE_UNCERTAINTY_UNKNOWN
```

A residual associated with stale or highly dispersed TLE references must receive:

```text
REFERENCE_LIMITED
```

and cannot reach the highest evidence levels.

---

# 20. Diagnostic Classes

Implement non-exclusive evidence tags:

```text
NORMAL
STATISTICAL_ANOMALY
STATION_SPECIFIC
CALIBRATION_CORRELATED
TIMING_CORRELATED
GEOMETRY_CORRELATED
LOCALISATION_CORRELATED
TRAJECTORY_MODEL_CORRELATED
REFERENCE_LIMITED
NETWORK_CONSISTENT
UNEXPLAINED
```

Multiple tags may coexist.

---

# 21. Evidence-Level Assignment

Suggested transparent rule-based logic:

```text
Level 0:
    no statistical anomaly

Level 1:
    statistical anomaly exists

Level 2:
    station-specific or strong leave-one-out influence

Level 3:
    instrument/calibration/timing/localisation correlation

Level 4:
    anomaly persists across station subsets
    AND not geometry-correlated
    AND not trajectory-model-correlated
    AND reference quality acceptable

Level 5:
    Level 4
    AND statistically strong
    AND persistent over multiple epochs
    AND no available diagnostic explanation
```

Do not use machine learning for the first implementation.

---

# 22. Statistical Anomaly Definition

Use Stage-6 uncertainty-aware statistics.

Preferred:

```text
position_mahalanobis_total_cov
```

when empirical TLE reference covariance exists.

Fallback:

```text
position_mahalanobis_obs_cov
```

with clear conditional interpretation.

Allow configurable:

```text
--anomaly-threshold
```

A default may correspond to a conservative chi-square threshold, but it must be documented and configurable.

---

# 23. Persistence Criterion

Single-epoch excursions are weaker evidence.

Define persistence metrics:

```text
consecutive anomalous epochs
anomalous fraction in rolling window
duration_s
```

Output:

```text
anomaly_run_length
anomaly_duration_s
```

Support:

```text
--min-persistent-epochs
```

for higher evidence levels.

---

# 24. Station Health Summary

Produce one summary per station.

Fields:

```text
station_id
n_epochs
usable_fraction
median_direction_residual_deg
median_line_residual_m
median_leave_one_out_shift_m
median_leave_one_out_mahalanobis
zero_gain_fraction
median_flagged_fraction
calibration_change_count
median_snr
median_image_visibility_offset_deg
timing_offset_median_s
timing_offset_mad_s
evidence_tags
health_flags
```

Do not collapse this into a single opaque health score by default.

---

# 25. Optional Transparent Health Score

If implemented, any health score must be transparent.

Possible components:

```text
direction residual percentile
line residual percentile
leave-one-out influence percentile
flagged fraction percentile
timing mismatch percentile
```

Output each component separately.

If a combined score is produced, document the formula exactly.

---

# 26. Epoch Diagnosis Table

Create:

```text
diagnosis.parquet
```

with one row per Stage-6 epoch.

Fields:

```text
time_utc
norad_id
position_residual_3d_m
dr_radial_m
dr_along_track_m
dr_cross_track_m
position_mahalanobis_obs_cov
position_mahalanobis_total_cov
stage5_nis
stage5_flags
stage4_condition_number
stage4_min_crossing_angle_deg
max_station_line_residual_m
max_station_direction_residual_deg
max_leave_one_out_shift_m
most_influential_station
tle_age_days
precision_validation
diagnostic_tags
evidence_level
diagnostic_confidence
quality_flags
```

`diagnostic_confidence` must be descriptive, not a fake calibrated probability unless explicitly calibrated.

---

# 27. Station Influence Table

Create:

```text
station_influence.parquet
```

Long format:

```text
time_utc
station_id
loo_x_ecef_m
loo_y_ecef_m
loo_z_ecef_m
loo_shift_3d_m
loo_shift_radial_m
loo_shift_along_track_m
loo_shift_cross_track_m
loo_mahalanobis
loo_condition_number
loo_min_crossing_angle_deg
degraded_geometry
direction_residual_angle_deg
direction_residual_az_deg
direction_residual_el_deg
line_residual_m
calibration_warning
timing_warning
stage3_flags
stage4_flags
```

---

# 28. Correlation Output

Write:

```text
correlations.csv
```

Fields:

```text
scope
station_id
x_variable
y_variable
method
coefficient
p_value
n_samples
interpretation_tag
```

Examples:

```text
zero_gain_fraction vs direction_residual_angle_deg
flagged_fraction vs line_residual_m
station_time_offset_s vs dr_along_track_m
condition_number vs position_residual_3d_m
NIS vs position_residual_3d_m
```

---

# 29. Stage-7 Summary

Write:

```text
summary.json
```

with:

```text
overall evidence level
overall diagnostic tags
number of anomalous epochs
persistent anomaly intervals
most influential station
station summaries
strongest calibration correlation
strongest timing correlation
strongest geometry correlation
strongest trajectory-model correlation
reference-quality status
precision-validation status
unexplained intervals
```

---

# 30. Stage-7 Manifest

Write:

```text
stage7_manifest.json
```

with:

```json
{
  "schema_version": 1,
  "stage": 7,
  "input_stage6_manifest": "...",
  "input_hashes": {},
  "alignment": {},
  "leave_one_out": {},
  "pairwise": {},
  "correlation_analysis": {},
  "evidence_rules": {},
  "quality_policy": {},
  "software": {},
  "provenance": {},
  "outputs": {}
}
```

Record exact rule thresholds.

---

# 31. Diagnostic Plots

Implement at minimum:

```text
RTN residuals with evidence levels
leave-one-out influence per station
station direction residual per station
strongest calibration correlations
timing vs along-track residual
geometry vs residual
station-by-time influence heatmap
```

Use matplotlib.

---

# 32. Quality Policy

Implement:

```text
strict
permissive
```

Default:

```text
strict
```

Strict policy should exclude rows already marked invalid by upstream stages.

Advisory flags must remain visible.

Do not manufacture diagnosis when the upstream chain is invalid.

---

# 33. Missing Data

Stage 7 must tolerate incomplete diagnostics.

For example, missing Stage-2 calibration metadata must not crash the analysis.

Instead record:

```text
CALIBRATION_METADATA_MISSING
```

and lower attribution confidence.

Do not treat missing evidence as evidence of absence.

---

# 34. Confidence Model

Implement a transparent rule-based confidence level using components such as:

```text
sample support
effect size
persistence
cross-method agreement
reference quality
data completeness
```

Map to:

```text
LOW
MODERATE
HIGH
```

Do not output numerical probabilities unless calibrated.

---

# 35. Synthetic Validation

Synthetic tests are mandatory.

## 35.1 Clean network

Inject consistent measurements across all stations.

Expected:

```text
Level 0
NORMAL
```

## 35.2 One biased station

Inject a fixed angular bias into one station.

Expected:

```text
STATION_SPECIFIC
large leave-one-out influence
Level >= 2
```

## 35.3 Calibration-correlated bias

Inject bias only when a synthetic calibration flag is active.

Expected:

```text
CALIBRATION_CORRELATED
Level >= 3
```

## 35.4 Timing-correlated bias

Inject timing offsets causing along-track displacement.

Expected:

```text
TIMING_CORRELATED
```

## 35.5 Poor geometry

Degrade crossing angles without instrument bias.

Expected:

```text
GEOMETRY_CORRELATED
```

not `NETWORK_CONSISTENT`.

## 35.6 Dynamics mismatch

Inject unmodelled acceleration.

Expected:

```text
TRAJECTORY_MODEL_CORRELATED
```

with Stage-5 innovation diagnostics.

## 35.7 Stale TLE

Use a deliberately poor reference.

Expected:

```text
REFERENCE_LIMITED
```

and cap maximum evidence level.

## 35.8 Network-consistent offset

Inject the same physically consistent trajectory offset seen by all stations.

Expected:

```text
NETWORK_CONSISTENT
```

with low station-specific influence.

## 35.9 Unexplained persistent residual

Create a persistent network-consistent residual with clean instrument, geometry, dynamics and reference diagnostics.

Expected:

```text
UNEXPLAINED
Level 5
```

---

# 36. Leave-One-Out Synthetic Test

For a known station bias, verify:

```text
full solution error is large
excluding biased station reduces error
excluding healthy stations does not
```

This test is mandatory.

---

# 37. Correlation False-Positive Control

Synthetic tests must include unrelated random metadata.

Verify that attribution logic does not systematically label random correlations as physical evidence.

Require:

```text
minimum sample count
effect-size threshold
multiple-testing caution in documentation
```

Do not promote diagnosis based on p-value alone.

---

# 38. Real-Data Acceptance Test

Use a complete Stage-6 run with:

```text
precision_validation = true
```

if available.

Run:

```bash
tart-diagnose \
    --stage6-manifest <stage6_manifest.json> \
    --output-dir results/stage7
```

Inspect:

```text
station influence
station direction biases
calibration correlations
timing correlations
geometry dependence
Stage-5 dynamics diagnostics
reference quality
evidence levels
```

If only stale-TLE Stage-6 data are available, preserve:

```text
REFERENCE_LIMITED
```

and do not classify the result as a high-confidence network anomaly.

---

# 39. Multi-Satellite Design

The Stage-7 code should be written so results from multiple satellites can later be aggregated.

Do not hard-code one NORAD ID.

The single-run output remains one satellite/window.

Future aggregation should be able to combine:

```text
satellite
epoch
station
diagnostic class
```

without schema changes.

---

# 40. No Automatic Physical Attribution

Do not create labels such as:

```text
SATELLITE_MANOEUVRE
CLOCK_FAILURE
RFI_EVENT
HARDWARE_FAILURE
```

unless a future stage adds independent evidence sufficient to support them.

Stage 7 should report evidence patterns, not unsupported physical causes.

---

# 41. API Design

Suggested interfaces:

```python
def load_diagnostic_chain(
    stage6_manifest: Path,
) -> DiagnosticInputs:
    """Load and verify linked outputs from Stages 2–6."""
```

```python
def compute_leave_one_out(
    inputs: DiagnosticInputs,
) -> pd.DataFrame:
    """Recompute localisation while excluding each station in turn."""
```

```python
def compute_station_direction_residuals(
    inputs: DiagnosticInputs,
) -> pd.DataFrame:
    """Compare each station direction with the network-derived Stage-5 trajectory."""
```

```python
def compute_correlations(
    table: pd.DataFrame,
    method: str = "spearman",
) -> pd.DataFrame:
    """Calculate diagnostic correlations with sample counts and significance metadata."""
```

```python
def assign_evidence(
    epoch: DiagnosticEpoch,
    rules: EvidenceRules,
) -> Diagnosis:
    """Assign diagnostic tags and evidence level using transparent rules."""
```

---

# 42. Logging

At INFO level report:

```text
Stage-6 input
NORAD ID
number of epochs
number of stations
number of valid aligned epochs
leave-one-out runs
pairwise runs
number of anomalous epochs
number of persistent intervals
most influential station
strongest diagnostic correlations
reference-quality status
maximum evidence level
output paths
```

Use `logging`.

---

# 43. Failure Behaviour

Add:

```text
Stage7Error
InvalidStage6ManifestError
ProvenanceChainError
AlignmentError
InfluenceAnalysisError
CorrelationAnalysisError
DiagnosisError
```

CLI exit codes:

```text
0 successful diagnosis
2 invalid input/configuration
3 completed but insufficient evidence/data for diagnosis
1 unexpected failure
```

---

# 44. Acceptance Criteria

## CLI

- [ ] `tart-diagnose` exists.
- [ ] Existing `tart-*` commands remain unchanged.
- [ ] `tests/test_cli_names.py` includes `tart-diagnose`.

## Functional

- [ ] Stage-6 output loads.
- [ ] Upstream provenance chain loads.
- [ ] Time alignment is explicit.
- [ ] Leave-one-station-out analysis works.
- [ ] Pairwise station diagnostics work.
- [ ] Station direction residuals are computed.
- [ ] Calibration metadata correlations are computed.
- [ ] Timing correlations are computed.
- [ ] Geometry dependence is computed.
- [ ] Stage-5 diagnostics are incorporated.
- [ ] TLE/reference quality is incorporated.
- [ ] Evidence tags are assigned.
- [ ] Evidence level is assigned.
- [ ] Station summaries are produced.
- [ ] Epoch diagnosis table is produced.
- [ ] Manifest and plots are produced.

## Scientific

- [ ] One biased station is identified in synthetic tests.
- [ ] Clean network remains Level 0.
- [ ] Poor geometry is not misclassified as network anomaly.
- [ ] Timing bias appears primarily along-track.
- [ ] Calibration-correlated synthetic bias is identified.
- [ ] Dynamics mismatch is identified through Stage-5 diagnostics.
- [ ] Stale TLE limits diagnosis confidence.
- [ ] Network-consistent synthetic offset is distinguished from station bias.
- [ ] Unexplained Level-5 case requires persistence and clean diagnostics.
- [ ] Missing metadata does not become false evidence.

## Engineering

- [ ] Stage-1 tests pass.
- [ ] Stage-2 tests pass.
- [ ] Stage-3 tests pass.
- [ ] Stage-4 tests pass.
- [ ] Stage-5 tests pass.
- [ ] Stage-6 tests pass.
- [ ] Stage-7 tests pass.
- [ ] Registry tests pass.
- [ ] No network access is required for unit tests.
- [ ] Public functions have docstrings.
- [ ] Logging is consistent.

---

# 45. Expected End State

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
tart-diagnose ...
```

Conceptually:

```text
TART visibilities
      |
direction measurements
      |
3-D localisation
      |
trajectory inference
      |
TLE comparison
      |
residuals
      |
tart-diagnose
      |
station influence
instrument correlations
timing correlations
geometry dependence
reference limitations
network-consistency evidence
```

---

# 46. Final Deliverables

The coding agent must deliver:

```text
src/tart_gnss/stage7/*
tests/stage7/*
docs/stage7_diagnosis.md
README.md update
pyproject.toml entry-point/dependency updates
tests/test_cli_names.py update
tart_gnss_stage7_implementation.md
```

Also provide:

```text
synthetic diagnosis example
real-data diagnosis example when suitable Stage-6 data exist
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

Run a synthetic end-to-end demonstration containing at least:

```text
clean network
one biased station
timing-correlated case
geometry-limited case
network-consistent unexplained case
```

Report:

```text
injected cause pattern
recovered diagnostic tags
evidence level
most influential station
leave-one-out shift
strongest correlation
whether the classification matched expectation
```

The implementation is complete only when Stage 7 can distinguish station/instrument/geometry/reference-correlated residuals from persistent network-consistent unexplained residuals without making unsupported physical-cause claims.
