# TART-GNSS Campaign Analysis Implementation Specification
## Multi-Satellite Population Statistics, Network Performance, and Scientific Reporting

**Repository:** `git@github.com:ifornax/tart-gnss.git`  
**Base branch:** `main`  
**Pipeline status before this work:** Stages 1–7 implemented  
**This specification implements:** a campaign/population-analysis layer across completed Stage-7 runs  
**Canonical CLI:** `tart-campaign`

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
tart-diagnose
```

This implementation must add:

```text
tart-campaign
```

Do not add public commands with any other naming pattern.

Update `pyproject.toml`, `README.md`, docs, `tests/test_cli_names.py`, and examples. Historical provenance files must remain unchanged.

---

# 1. Purpose

The single-pass pipeline is complete through Stage 7.

The scientific question now changes from:

```text
Can the TART network reconstruct and diagnose one satellite pass?
```

to:

```text
How well does the TART network perform statistically across many satellites,
epochs, stations, geometries, constellations, and observing conditions?
```

The campaign layer must aggregate multiple completed Stage-7 runs and produce:

1. population-level localisation and trajectory statistics;
2. network performance versus observing conditions;
3. station stability and station-to-station comparisons;
4. constellation and satellite comparisons;
5. residual distributions in Cartesian and RTN coordinates;
6. anomaly/evidence-level occurrence rates;
7. campaign-level quality and reproducibility summaries;
8. scientifically defensible tables and plots suitable for a talk or paper.

This layer is **not** another single-pass processing stage.

It must consume completed pipeline products rather than recomputing Stages 1–7.

---

# 2. Scientific Questions

The implementation must support questions such as:

```text
What is the median angular localisation accuracy across the campaign?

What is the median 3-D localisation uncertainty?

How do radial, along-track, and cross-track residuals vary across satellites?

How does performance depend on:
    elevation,
    SNR,
    station geometry,
    crossing angle,
    condition number,
    flagged fraction,
    calibration state,
    timing offsets?

Which stations are most stable?

Are some stations systematically biased?

Do residuals differ between GPS, Galileo, BeiDou, etc.?

How often are anomalies:
    station-specific,
    calibration-correlated,
    timing-correlated,
    geometry-correlated,
    trajectory-model-correlated,
    reference-limited,
    network-consistent,
    unexplained?

What fraction of passes are suitable for precision validation?

How reproducible are results across repeated passes of the same satellite?
```

---

# 3. Input Philosophy

The campaign layer must operate from completed Stage-7 outputs.

Primary per-run input:

```text
stage7_manifest.json
```

The campaign tool must follow each Stage-7 provenance chain backward only as needed.

Do not rerun Stage 1–7.

Do not recompute scientific products already present upstream unless required for a campaign-derived statistic that cannot otherwise be calculated.

---

# 4. Input Discovery

Support two modes.

## 4.1 Explicit manifest list

```bash
tart-campaign \
    --stage7-manifest run1/stage7_manifest.json \
    --stage7-manifest run2/stage7_manifest.json \
    --stage7-manifest run3/stage7_manifest.json
```

## 4.2 Campaign root discovery

```bash
tart-campaign \
    --input-root results/stage7
```

When using `--input-root`, recursively discover only:

```text
stage7_manifest.json
```

Do not infer runs from arbitrary CSV or Parquet files.

Every discovered run must be validated through its manifest.

---

# 5. CLI

Implement:

```bash
tart-campaign
```

Example:

```bash
tart-campaign \
    --input-root results/stage7 \
    --output-dir results/campaign \
    --campaign-name gnss_validation_2026
```

Required/supported options:

```text
--input-root
--stage7-manifest
--output-dir
--campaign-name
--quality-policy
--precision-only
--min-epochs-per-run
--min-runs-per-group
--group-by
--reference-frame
--bootstrap
--seed
--no-plots
--verbose
```

Useful optional filters:

```text
--norad-id
--constellation
--station
--start
--end
--min-elevation
--max-tle-age-days
--max-evidence-level
--precision-validation-only
```

At least one of `--input-root` or one-or-more `--stage7-manifest` inputs is required.

---

# 6. Canonical Campaign Unit

Define one campaign observation unit as:

```text
one completed Stage-7 run
=
one satellite × one common observing window
```

Each run must have a unique campaign identifier.

Recommended:

```text
run_id
```

constructed from date, NORAD ID, and window index/start time, while preserving original manifest identity where available.

---

# 7. Repository Integration

Add:

```text
src/tart_gnss/campaign/
├── __init__.py
├── cli.py
├── discovery.py
├── loader.py
├── schema.py
├── filters.py
├── aggregation.py
├── accuracy.py
├── station_metrics.py
├── constellation_metrics.py
├── anomaly_metrics.py
├── reproducibility.py
├── statistics.py
├── bootstrap.py
├── output.py
├── plotting.py
└── errors.py
```

Tests:

```text
tests/campaign/
├── __init__.py
├── test_discovery.py
├── test_loader.py
├── test_schema.py
├── test_filters.py
├── test_aggregation.py
├── test_accuracy.py
├── test_station_metrics.py
├── test_constellation_metrics.py
├── test_anomaly_metrics.py
├── test_reproducibility.py
├── test_statistics.py
├── test_bootstrap.py
├── test_output.py
├── test_cli.py
└── test_synthetic_campaign.py
```

Documentation:

```text
docs/campaign_analysis.md
```

---

# 8. Campaign Provenance

Every loaded Stage-7 manifest must record:

```text
path
SHA-256
git commit
pipeline version
NORAD ID
observation date
time span
station set
TLE source
TLE epoch
precision-validation status
overall evidence level
```

Campaign output must retain the exact list of included and excluded runs.

Do not silently discard invalid runs.

---

# 9. Schema Harmonisation

The loader must support the current Stage-7 schema version and fail clearly on incompatible future versions.

Create a canonical internal campaign schema.

At minimum, per run:

```text
run_id
norad_id
satellite_name
constellation
start_utc
end_utc
duration_s
n_epochs
n_valid_epochs
stations
n_stations
precision_validation
tle_epoch_utc
median_tle_age_days
max_abs_tle_age_days
overall_evidence_level
overall_diagnostic_tags
diagnostic_confidence
```

---

# 10. Per-Epoch Campaign Table

Create:

```text
campaign_epochs.parquet
```

Recommended columns:

```text
run_id
time_utc
norad_id
constellation
position_residual_3d_m
velocity_residual_3d_m_s
dr_radial_m
dr_along_track_m
dr_cross_track_m
position_mahalanobis_obs_cov
position_mahalanobis_total_cov
stage4_sigma_3d_m
stage5_sigma_3d_m
stage4_condition_number
stage4_min_crossing_angle_deg
stage4_max_crossing_angle_deg
stage5_nis
evidence_level
diagnostic_tags
precision_validation
tle_age_days
```

Add upstream station/network quality features when available.

---

# 11. Per-Run Summary Table

Create:

```text
campaign_runs.parquet
```

One row per Stage-7 run.

Fields should include:

```text
run_id
norad_id
satellite_name
constellation
start_utc
duration_s
stations
n_stations
n_epochs
n_valid_epochs
usable_fraction
median_elevation_deg
median_snr
median_flagged_fraction
median_condition_number
median_min_crossing_angle_deg
median_stage4_sigma_3d_m
median_stage5_sigma_3d_m
median_position_residual_3d_m
rms_position_residual_3d_m
median_radial_m
median_along_track_m
median_cross_track_m
rms_radial_m
rms_along_track_m
rms_cross_track_m
median_velocity_residual_m_s
median_position_mahalanobis_obs_cov
median_position_mahalanobis_total_cov
precision_validation
median_tle_age_days
overall_evidence_level
overall_diagnostic_tags
diagnostic_confidence
most_influential_station
```

---

# 12. Station-Level Campaign Table

Create:

```text
campaign_stations.parquet
```

One row per:

```text
run_id × station
```

Fields:

```text
run_id
station_id
norad_id
constellation
n_epochs
usable_fraction
median_direction_residual_deg
median_line_residual_m
median_leave_one_out_shift_m
median_leave_one_out_mahalanobis
median_snr
median_flagged_fraction
zero_gain_fraction
calibration_change_count
median_time_offset_s
mad_time_offset_s
station_tags
health_flags
```

---

# 13. Satellite Metadata

Determine constellation from existing pipeline metadata.

Supported canonical values:

```text
GPS
Galileo
BeiDou
GLONASS
Unknown
```

Do not infer constellation from satellite name with fragile string parsing when structured metadata exists.

Preserve NORAD ID, satellite name, and constellation.

---

# 14. Accuracy Terminology

Be precise with terminology.

The campaign must distinguish:

```text
measurement uncertainty
reference residual
repeatability
accuracy relative to TLE
```

Do not call the median observed-minus-TLE residual an absolute astrometric accuracy unless reference quality supports that interpretation.

Preferred wording:

```text
TLE-relative residual
trajectory residual
network repeatability
measurement uncertainty
```

---

# 15. Angular Localisation Performance

Aggregate Stage-3 directional performance.

Per station and campaign:

```text
median sigma_az_deg
median sigma_el_deg
median angular covariance scale
median visibility-image disagreement_deg
```

If an observed-vs-reference direction residual is available, keep it distinct from formal fit uncertainty.

---

# 16. 3-D Localisation Performance

Aggregate Stage-4 metrics:

```text
median 3-D uncertainty
median sigma_x
median sigma_y
median sigma_z
covariance eigenvalues
major-axis uncertainty
minor-axis uncertainty
```

Prefer covariance-eigenvalue summaries because the uncertainty is strongly anisotropic.

Output:

```text
median_sigma_major_m
median_sigma_intermediate_m
median_sigma_minor_m
```

---

# 17. Trajectory Performance

Aggregate Stage-5 metrics:

```text
filtered uncertainty
smoothed uncertainty
position RMSE proxy relative to Stage-4
velocity uncertainty
median NIS
NIS exceedance fraction
innovation autocorrelation
filter reset count
dynamics mismatch fraction
```

Do not treat TLE residual as Stage-5 accuracy.

---

# 18. TLE-Relative Residual Statistics

For runs with:

```text
precision_validation = true
```

aggregate:

```text
3-D residual norm
radial residual
along-track residual
cross-track residual
velocity residual
```

Report:

```text
median
mean
RMS
MAD
16th percentile
84th percentile
95th percentile
```

Do not aggregate stale-TLE runs into precision-performance metrics by default.

---

# 19. Precision and Descriptive Cohorts

Every campaign run must belong to one of:

```text
precision_validation
descriptive_only
```

Default scientific summaries must use `precision_validation = true` unless explicitly overridden.

Still report counts of excluded descriptive runs.

---

# 20. TLE Age Dependence

Analyse residuals versus absolute and signed TLE age.

Suggested bins:

```text
0–1 day
1–3 days
3–7 days
7–14 days
>14 days
```

Report `n`, median residual, RMS residual, and MAD per bin.

Do not infer causation automatically.

---

# 21. Elevation Dependence

Suggested bins:

```text
20–30 deg
30–40 deg
40–50 deg
50–60 deg
60–70 deg
70–90 deg
```

Metrics:

```text
direction uncertainty
direction residual
Stage-4 uncertainty
Stage-6 residual
```

Use station-level elevation where appropriate.

---

# 22. SNR Dependence

Aggregate:

```text
direction error/uncertainty vs SNR
Stage-4 uncertainty vs SNR
station influence vs SNR
```

Use logarithmic SNR binning only if justified by the data.

---

# 23. Geometry Dependence

Analyse performance versus:

```text
minimum crossing angle
condition number
number of stations
```

At minimum produce:

```text
median Stage-4 3-D uncertainty vs crossing angle
median residual vs crossing angle
median leave-one-out shift vs crossing angle
```

---

# 24. Station Stability

Across all runs, calculate for each station:

```text
n_runs
n_epochs
usable_fraction
median direction residual
MAD direction residual
median line residual
median leave-one-out influence
median flagged fraction
zero-gain incidence
gain-change incidence
timing-offset median
timing-offset MAD
fraction of runs tagged:
    STATION_SPECIFIC
    CALIBRATION_CORRELATED
    TIMING_CORRELATED
    LOCALISATION_CORRELATED
```

This table is a key campaign output.

---

# 25. Station Ranking

Do not create a simplistic “best station” score by default.

Instead produce a transparent station comparison table.

If a ranking is implemented, expose each component and all weights.

---

# 26. Constellation Comparison

Group results by:

```text
GPS
Galileo
BeiDou
GLONASS
Unknown
```

Compare:

```text
number of runs
median SNR
median direction uncertainty
median Stage-4 uncertainty
median TLE-relative residual
evidence-level distribution
precision-validation fraction
```

Default minimum group size:

```text
--min-runs-per-group 3
```

Flag smaller groups as under-sampled.

---

# 27. Satellite-Level Comparison

For each NORAD ID:

```text
n_runs
time span
median residual
run-to-run scatter
median uncertainty
evidence-level distribution
```

This supports repeatability analysis.

---

# 28. Repeatability Across Passes

For satellites observed multiple times, calculate per-pass RTN summaries and then:

```text
between-pass standard deviation
between-pass MAD
within-pass median uncertainty
```

Compare between-pass scatter with within-pass uncertainty.

---

# 29. Run-to-Run Reproducibility

For repeated observations with similar geometry, quantify:

```text
direction-bias repeatability
station influence repeatability
RTN residual repeatability
evidence-tag repeatability
```

Aggregate per pass rather than per epoch.

---

# 30. Evidence-Level Population Statistics

Calculate separately by run, epoch, and station:

```text
fraction Level 0
fraction Level 1
fraction Level 2
fraction Level 3
fraction Level 4
fraction Level 5
```

Keep denominators separate.

---

# 31. Diagnostic Tag Frequencies

Calculate frequencies for:

```text
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

Output count, fraction, and confidence interval where appropriate.

---

# 32. Confidence Intervals

For medians and other non-Gaussian population metrics, use bootstrap confidence intervals.

Support:

```text
--bootstrap N
```

Suggested default:

```text
1000
```

for campaign-level reports, while remaining configurable.

Bootstrap at the **run level**, not the epoch level.

---

# 33. Cluster-Aware Statistics

Epochs within one pass are correlated.

Therefore top-level campaign inference must use:

```text
run as the independent cluster
```

For grouped comparisons, aggregate per run first.

Do not report tiny p-values caused by treating thousands of correlated epochs as independent observations.

---

# 34. Statistical Comparisons

For two groups, support:

```text
Mann-Whitney U
bootstrap difference in medians
```

For more than two groups:

```text
Kruskal-Wallis
```

Use only when group sample sizes are adequate.

Report effect sizes alongside p-values.

---

# 35. Regression / Dependence Models

Optional campaign-level regressions may relate performance to:

```text
elevation
SNR
crossing angle
condition number
TLE age
```

Recommended initial approaches:

```text
Theil-Sen
robust linear regression
Spearman correlation
```

Do not introduce complex machine learning in the initial implementation.

---

# 36. Anomaly Rate versus Conditions

Calculate diagnostic-tag rates as functions of station, constellation, elevation, geometry, calibration state, and TLE age.

Examples:

```text
fraction STATION_SPECIFIC vs station
fraction GEOMETRY_CORRELATED vs crossing-angle bin
fraction REFERENCE_LIMITED vs TLE-age bin
```

---

# 37. Campaign Quality Filters

Implement:

```text
strict
permissive
```

Default:

```text
strict
```

Strict mode should exclude runs with broken provenance or insufficient valid epochs from scientific summaries.

Still list excluded runs separately.

---

# 38. Exclusion Reasons

Create:

```text
excluded_runs.csv
```

Fields:

```text
run_id
manifest_path
norad_id
reason
quality_flags
```

Example reasons:

```text
BROKEN_PROVENANCE
INSUFFICIENT_EPOCHS
NO_PRECISION_VALIDATION
MISSING_STAGE7_OUTPUT
SCHEMA_MISMATCH
```

---

# 39. Campaign Output Layout

Recommended:

```text
results/campaign/<campaign_name>/
├── campaign_manifest.json
├── campaign_runs.parquet
├── campaign_epochs.parquet
├── campaign_stations.parquet
├── excluded_runs.csv
├── summary.json
├── tables/
│   ├── station_summary.csv
│   ├── constellation_summary.csv
│   ├── satellite_summary.csv
│   ├── evidence_summary.csv
│   ├── tle_age_summary.csv
│   └── geometry_summary.csv
├── diagnostics/
│   ├── residual_distribution.png
│   ├── rtn_distribution.png
│   ├── residual_vs_elevation.png
│   ├── residual_vs_snr.png
│   ├── uncertainty_vs_crossing_angle.png
│   ├── station_stability.png
│   ├── constellation_comparison.png
│   ├── evidence_level_distribution.png
│   ├── tle_age_dependence.png
│   └── repeatability.png
└── logs/
    └── campaign.log
```

---

# 40. Campaign Manifest

Write:

```text
campaign_manifest.json
```

with:

```json
{
  "schema_version": 1,
  "analysis": "campaign",
  "campaign_name": "...",
  "input_manifests": [],
  "input_hashes": {},
  "filters": {},
  "grouping": {},
  "bootstrap": {},
  "quality_policy": {},
  "software": {},
  "provenance": {},
  "outputs": {}
}
```

Record all included/excluded runs, manifest hashes, CLI command, git commit, software versions, random seed, bootstrap count, and filter settings.

---

# 41. Summary JSON

Write:

```text
summary.json
```

At minimum:

```text
campaign_name
n_discovered_runs
n_included_runs
n_excluded_runs
n_satellites
n_constellations
n_unique_stations
observation_start
observation_end
total_observing_duration_s
precision_validation_fraction
median_stage3_direction_uncertainty
median_stage4_3d_uncertainty
median_stage5_smoothed_uncertainty
median_position_residual_3d
median_radial_residual
median_along_track_residual
median_cross_track_residual
bootstrap confidence intervals
evidence-level fractions
diagnostic-tag fractions
most stable stations
most frequently flagged stations
repeatability summary
```

Use “TLE-relative residual” rather than absolute accuracy unless justified.

---

# 42. Core Scientific Plots

Implement at minimum:

1. 3-D residual distribution / ECDF.
2. RTN residual distributions.
3. Residual versus elevation.
4. Residual versus SNR.
5. Stage-4 uncertainty versus crossing angle.
6. Station stability multi-panel summary.
7. Constellation comparison when adequately sampled.
8. Evidence-level distribution.
9. TLE-age dependence.
10. Repeatability across passes.

Use matplotlib only.

---

# 43. Presentation-Ready Figures

Add:

```text
--publication-plots
```

When enabled, produce larger fonts, vector PDF + PNG, descriptive axis labels, units, and talk/paper-friendly legends.

Do not hard-code stylistic colours unless an existing project convention specifies them.

---

# 44. Scientific Reporting Tables

Generate:

```text
station_summary.csv
constellation_summary.csv
satellite_summary.csv
evidence_summary.csv
geometry_summary.csv
tle_age_summary.csv
```

Each table should contain, where relevant:

```text
n
median
MAD
16th percentile
84th percentile
bootstrap CI
```

---

# 45. Campaign Comparison Mode

An optional future extension may compare two campaign manifests.

If implemented, compare residual distributions, station stability, evidence-level frequencies, and precision-validation fraction.

This feature may be deferred if it complicates the first implementation.

---

# 46. No Reprocessing Rule

`tart-campaign` must not invoke:

```text
tart-localise
tart-triangulate
tart-track
tart-validate
tart-diagnose
```

It consumes completed outputs only.

If a run is incomplete, report it as incomplete.

Do not silently finish the pipeline.

---

# 47. Synthetic Campaign Validation

Synthetic tests are mandatory.

Create multiple synthetic Stage-7 runs with known population properties.

## 47.1 Clean campaign

Expected:

```text
low residuals
mostly Level 0
stable stations
```

## 47.2 One unstable station

Inject repeated station bias into one station across multiple runs.

Expected:

```text
higher station residual
higher STATION_SPECIFIC fraction
larger LOO influence
```

## 47.3 Geometry-dependent campaign

Inject varying crossing angles.

Expected:

```text
Stage-4 uncertainty worsens as crossing angle decreases
```

## 47.4 TLE-age campaign

Inject residual growth with TLE age.

Expected positive dependence recovered.

## 47.5 Constellation difference

Inject different SNR/noise distributions by constellation and recover the grouped difference.

## 47.6 Repeatability

Repeated passes of the same synthetic satellite should recover known pass-to-pass scatter.

---

# 48. Bootstrap Validation

Verify run-level bootstrap confidence intervals using synthetic known populations.

Test coverage of median residual, coverage of station-summary statistics, and reproducibility with fixed random seed.

---

# 49. Correlation Safety

Campaign-level correlations must not use all epochs as independent samples when evaluating significance.

Use per-run aggregates or cluster-aware/bootstrap approaches.

Document this explicitly.

---

# 50. Real Campaign Acceptance Test

Use all completed real Stage-7 runs available locally.

Run:

```bash
tart-campaign \
    --input-root results/stage7 \
    --campaign-name real_tart_gnss \
    --output-dir results/campaign
```

Report:

```text
number discovered
number included
number precision-valid
satellites
stations
constellations
time span
```

Do not claim scientific population results if only one or two runs exist.

The tool should still generate descriptive outputs while flagging:

```text
INSUFFICIENT_CAMPAIGN_SIZE
```

---

# 51. Campaign Quality Flags

Suggested flags:

```text
GOOD
INSUFFICIENT_CAMPAIGN_SIZE
INSUFFICIENT_GROUP_SIZE
LOW_PRECISION_VALIDATION_FRACTION
BROKEN_PROVENANCE_PRESENT
DOMINATED_BY_ONE_SATELLITE
DOMINATED_BY_ONE_STATION_SET
TLE_AGE_IMBALANCE
MISSING_CONSTELLATION_METADATA
```

---

# 52. Minimum Scientific Sample Size

Use transparent warnings rather than one hard law.

Suggested guidance:

```text
< 3 runs: descriptive only
3–9 runs: preliminary population statistics
>= 10 runs: basic campaign inference permitted
>= 30 runs: stronger grouped comparisons possible
```

---

# 53. Campaign Balance Diagnostics

Report whether the dataset is dominated by one satellite, one constellation, one station combination, one day, or one elevation range.

---

# 54. Station-Combination Analysis

Group runs by station set.

Compare:

```text
geometry
uncertainty
residual
evidence-level distribution
```

This can inform future network design.

---

# 55. Network Design Insight

Provide campaign summaries useful for deciding which station combinations are geometrically strongest.

Metrics:

```text
median minimum crossing angle
median condition number
median Stage-4 uncertainty
usable fraction
```

Do not automatically recommend removing a station solely because it performs poorly in one campaign.

---

# 56. No Automatic Physical Claims

Campaign analysis must not turn evidence tags into physical-cause claims.

Do not produce:

```text
satellite manoeuvre rate
clock failure rate
hardware failure rate
RFI event rate
```

unless independent evidence exists.

Report network-consistent unexplained fraction, station-specific fraction, and reference-limited fraction instead.

---

# 57. API Design

Suggested interfaces:

```python
def discover_stage7_runs(input_root: Path) -> list[Path]:
    """Discover Stage-7 manifests recursively."""
```

```python
def load_campaign_run(manifest_path: Path) -> CampaignRun:
    """Load and validate one completed Stage-7 run."""
```

```python
def build_campaign_tables(runs: list[CampaignRun]) -> CampaignTables:
    """Build run, epoch and station population tables."""
```

```python
def summarize_stations(station_table: pd.DataFrame) -> pd.DataFrame:
    """Compute cross-run station stability statistics."""
```

```python
def bootstrap_run_level(
    run_table: pd.DataFrame,
    statistic: Callable,
    n_bootstrap: int,
    seed: int,
) -> BootstrapResult:
    """Bootstrap campaign statistics using runs as the resampling unit."""
```

---

# 58. Logging

At INFO level report:

```text
campaign name
input root
manifests discovered
manifests loaded
runs included
runs excluded
precision-valid runs
unique satellites
constellations
stations
observation time span
bootstrap count
main summary metrics
output paths
```

Use `logging`.

---

# 59. Failure Behaviour

Add:

```text
CampaignError
ManifestDiscoveryError
InvalidStage7ManifestError
CampaignSchemaError
InsufficientCampaignDataError
AggregationError
BootstrapError
```

CLI exit codes:

```text
0  campaign analysis completed
2  invalid input/configuration
3  completed but insufficient data for scientific campaign inference
1  unexpected failure
```

The tool may still write descriptive outputs with exit code 3.

---

# 60. Acceptance Criteria

## CLI

- [ ] `tart-campaign` exists.
- [ ] Existing `tart-*` commands remain unchanged.
- [ ] `tests/test_cli_names.py` includes `tart-campaign`.

## Functional

- [ ] Explicit Stage-7 manifest input works.
- [ ] Recursive input-root discovery works.
- [ ] Stage-7 manifest validation works.
- [ ] Run-level table is produced.
- [ ] Epoch-level table is produced.
- [ ] Station-level table is produced.
- [ ] Excluded-runs table is produced.
- [ ] Precision/descriptive cohorts are separated.
- [ ] Station summaries are produced.
- [ ] Satellite summaries are produced.
- [ ] Constellation summaries are produced.
- [ ] Geometry summaries are produced.
- [ ] TLE-age summaries are produced.
- [ ] Evidence-level summaries are produced.
- [ ] Repeatability metrics are produced.
- [ ] Bootstrap CIs are produced.
- [ ] Campaign manifest is written.
- [ ] Diagnostic plots are produced.

## Scientific

- [ ] Top-level inference treats runs, not epochs, as independent units.
- [ ] Stale/reference-limited runs are excluded from precision metrics by default.
- [ ] Residual is not mislabeled as absolute accuracy.
- [ ] Station instability is recoverable in synthetic campaigns.
- [ ] Geometry dependence is recoverable.
- [ ] TLE-age dependence is recoverable.
- [ ] Repeated-pass scatter is correctly estimated.
- [ ] Group comparisons enforce minimum sample size.
- [ ] Effect sizes are reported with p-values.
- [ ] Campaign balance diagnostics are reported.

## Engineering

- [ ] Stages 1–7 tests pass.
- [ ] Campaign tests pass.
- [ ] Registry tests pass.
- [ ] No network access required for unit tests.
- [ ] Public functions have docstrings.
- [ ] Logging is consistent.
- [ ] Large observational data are not committed.

---

# 61. Expected End State

The complete workflow becomes:

```text
Per-pass pipeline:
    tart-common-gnss
    tart-retrieve
    tart-localise
    tart-triangulate
    tart-track
    tart-validate
    tart-diagnose

Campaign layer:
    tart-campaign
```

Conceptually:

```text
many completed Stage-7 runs
        |
        v
tart-campaign
        |
        +--> localisation performance
        +--> trajectory/TLE residual population
        +--> station stability
        +--> constellation comparison
        +--> geometry dependence
        +--> TLE-age dependence
        +--> evidence-level statistics
        +--> repeatability
        |
        v
scientific tables + presentation/paper figures
```

---

# 62. Final Deliverables

The coding agent must deliver:

```text
src/tart_gnss/campaign/*
tests/campaign/*
docs/campaign_analysis.md
README.md update
pyproject.toml entry-point/dependency updates
tests/test_cli_names.py update
tart_gnss_campaign_analysis_implementation.md
```

Also provide:

```text
synthetic multi-run campaign example
real campaign example using available Stage-7 runs
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

Run a synthetic campaign demonstration with at least:

```text
multiple satellites
multiple passes
multiple stations
varying geometry
one unstable station
multiple constellations
varying TLE age
```

Report:

```text
number of runs
number of satellites
precision-valid fraction
median 3-D localisation uncertainty
median TLE-relative residual
median radial/along-track/cross-track residual
station stability summary
constellation summary
geometry dependence
TLE-age dependence
evidence-level fractions
repeatability statistics
bootstrap confidence intervals
```

The implementation is complete only when the campaign layer can turn many independent Stage-7 runs into statistically defensible population-level network-performance results without treating correlated epochs as independent observations or overstating TLE-relative residuals as absolute truth.
