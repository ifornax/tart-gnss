# Stage 7: sensor-network diagnosis and systematic-effect evidence (`tart-diagnose`)

Stages 3–6 work out where the satellite is, how uncertain that estimate is,
and how the measured trajectory differs from a TLE. Stage 7 asks what
evidence says about **where a discrepancy comes from**:

- one station;
- an instrument or calibration effect;
- timing;
- network geometry;
- the trajectory model;
- the TLE reference;
- or something the network sees consistently and nothing explains.

It reports *evidence patterns*. It never asserts a physical cause: there are
no labels such as `SATELLITE_MANOEUVRE`, `CLOCK_FAILURE`, `RFI_EVENT` or
`HARDWARE_FAILURE`, and a test enforces this. It never modifies an upstream
output, and a test hashes every upstream file before and after.

Specification: [`tart_gnss_stage7_implementation.md`](../tart_gnss_stage7_implementation.md).

## 1. Inputs and provenance

The single input is `stage6_manifest.json`. Stage 7 follows the paths
**recorded** in each manifest:

```text
Stage 6 manifest ──input_stage5_manifest──▶ Stage 5 ──input_stage4_manifest──▶ Stage 4
        ──input_stage3_manifest──▶ Stage 3 ──input_manifest──▶ Stage 2
```

Each recorded path is resolved in this order:

1. as recorded;
2. relative to that run's `provenance.working_directory`;
3. next to the downstream manifest.

Where a downstream manifest recorded the upstream SHA-256, it is verified.
Every link is reported as `verified`, `unverified`, `mismatch` or `missing`.
Unrelated files are never searched for.

| Situation | Effect |
|---|---|
| hash mismatch under `strict` | fatal (exit 2) |
| hash mismatch under `permissive` | recorded (`PROVENANCE_UNVERIFIED`), data used |
| Stage 5 or Stage 4 missing | the chain is broken (exit 2) |
| Stage 3 missing | no re-triangulation |
| Stage 2 missing | `CALIBRATION_METADATA_MISSING` (§9) |

Loaded per stage:

| Stage | Tables |
|---|---|
| 6 | `residuals`, `tle_trajectory`, `summary.json` |
| 5 | `trajectory`, `innovations`, manifest run flags and resets |
| 4 | `positions`, `station_contributions`, frozen station coordinates |
| 3 | `directions` (SNR, reduced χ², flagged fraction, visibility–image separation, baselines) |
| 2 | per native file: time span, zero-gain antennas, gain-solution fingerprint |

Each station sample is given the Stage-2 calibration of the native file that
covers its time. That yields a zero-gain count and a gain-solution index,
which increases by one at each fingerprint change.

## 2. Alignment

Nothing is joined by row number. Each join uses the nearest time within
`--alignment-tolerance-s` (default 0.01 s):

- **Epoch table:** Stage-6 row ↔ Stage-5 trajectory and innovation row ↔
  Stage-4 position row ↔ TLE state, by `time_utc`.
- **Station table:** each Stage-4 station contribution ↔ its Stage-3 row, by
  (`station_id`, `station_time_utc`) ↔ the covering Stage-2 file.

Missing links are counted in the manifest. A test shuffles the upstream rows
and gets an identical alignment.

## 3. Leave-one-station-out and pairwise analysis

All solutions reuse the Stage-4 code:

- sight lines come from the Stage-3 rows via `stage4.loader.sight_line`, with
  the frozen coordinates;
- lines are intersected with `stage4.triangulation.triangulate_weighted`.

No triangulation maths is duplicated. As a check, the full solution is first
recomputed. On every synthetic chain and on the real data it reproduces the
Stage-4 table **exactly (0 m)**.

**Leave-one-out**, at each epoch with ≥ 3 stations and for each station i:

- `shift = x_full − x_−i`, in ECEF and in RTN (on the Stage-6 TLE basis);
- `loo_mahalanobis = shiftᵀ (C_−i − C_full)⁺ shift`.

For nested least-squares estimates, the covariance of their difference is
exactly C_−i − C_full. One station carries only two angular constraints, so
that matrix has rank ≤ 2. Its pseudo-inverse is used, and `loo_dof` records
the rank: the statistic is χ²_dof under the null. If the matrix is not
positive semidefinite, only the physical shift is reported.

With three stations, removing one leaves two. Those rows have
`degraded_geometry = true`, and their LOO condition number and crossing
angle are reported.

**Three-station caveat.** Removing a *healthy* station keeps a biased line in
the solution, so every LOO shift can be large. In the biased-station demo
the median shift is 39.5 / 104 / 229 km when removing Madagascar / Mauritius
/ South Africa. LOO influence alone therefore cannot name a culprit, and
neither can `most_influential_station`. That field is filled only when ≥ 3
epochs have a LOO analysis, and it is not evidence by itself.

**Pairwise**, at each epoch and for each station pair, the analysis computes:

- the two-line solution (`pair_position`);
- the closest-approach separation `pair_line_separation_m`;
- the line angle;
- the normalised separation
  `pair_separation_z = sep / √(nᵀ(ρ₁²C₁ + ρ₂²C₂)n)`, where n is the common
  normal and C the tangent angular covariances.

## 4. Station direction residuals

The Stage-5 smoothed state is moved to the station's own sample time
(`r + v dt`), rotated to ECEF and seen from the station. It is then compared
with the measured direction:

- the angle between them;
- (l, m) direction-cosine differences in the station's East/North frame,
  which stay defined at the zenith, unlike azimuth;
- Az/El differences;
- the components along and across the predicted apparent motion.

All are measured − predicted. Per station, the robust bias is the median and
MAD of the Az, El, angle and along-motion residuals, with descriptive slopes
against elevation and time. Biases are never calibrated away.

Because the comparison uses each station's own, correctly recorded time, a
station that samples late shows **no** direction residual. The effect
appears in the Stage-4 solution, which matches stations without
interpolation: in the line residuals and in the Stage-4 − TLE position. See
§5 (timing).

## 5. Correlation analysis

Each test pairs a variable from one *family* with a residual measure.

**Per station, over epochs** (`scope = station`):

| Family | Variables |
|---|---|
| calibration | zero-gain count, gain-solution index, flagged fraction, baselines used |
| localisation | SNR, reduced χ², visibility–image separation |
| timing | station time offset, and its absolute value |

The residual measures are the direction residual (angle and along-motion),
the line residual and the LOO shift.

**Per epoch** (`scope = network`):

| Family | Variables |
|---|---|
| geometry | predicted condition number, 1/sin(minimum crossing angle) |
| trajectory model | Stage-5 NIS |
| timing | largest and mean station time offset |

The residual measures are the Stage-4 − TLE distance, the Stage-6 3-D and
|along-track| residuals, and the largest LOO shift.

**Guarding against spurious "causes".** Within one pass, geometry, residuals
and metadata all drift smoothly with time, and two monotone series correlate
almost perfectly with no causal link. The method therefore:

1. computes the Spearman coefficient (default; `pearson` and `kendall` are
   available) and also a **partial** rank correlation after removing a
   linear trend in time. Decisions use the partial coefficient.
2. bases p-values on an effective sample size
   n(1 − r_x r_y)/(1 + r_x r_y), where r is the lag-1 autocorrelation of each
   detrended series;
3. requires `n_samples ≥ --min-samples` (8);
4. requires variation: ≥ 3 distinct values, or a two-valued indicator with ≥
   3 samples per group. Geometry variables must also vary over ≥ 5 % of
   their median;
5. requires an effect size `|partial| ≥ --min-correlation` (0.5);
6. applies a Holm-adjusted p < 0.01 within each family. A family holds about
   ten tests, so a single test needs p below roughly 10⁻³.

P-values are descriptive. A run makes about 130 tests, so some correlations
will look strong by chance (multiple testing). A p-value alone never
promotes a diagnosis.

Geometry variables are computed from the **predicted** sight lines (stations
to the Stage-5 position), not the measured ones. The measured condition
number carries the measurement noise. In the biased-station scenario, that
noise produced a spurious 0.93 "geometry correlation".

"Poorer geometry → larger Stage-4 σ" is labelled `EXPECTED`, never a tag.

**False-positive control (§37).** In 20 clean networks with random-walk
flagged-fraction and SNR metadata, no correlation is supported and no
evidence level is above 0. Ten seeds of random-walk metadata in the unit
tests give the same result.

## 6. Evidence

**Statistical anomaly.** The statistic is Stage-6
`position_mahalanobis_total_cov` when an empirical TLE covariance exists,
otherwise `position_mahalanobis_obs_cov` (conditional on the TLE). It is
compared with a threshold of χ²₃ at 1 − 0.0027/n, a Bonferroni correction
over the n diagnosed epochs. That gives 21.7 for 36 epochs, and
`--anomaly-threshold` overrides it.

**Per-epoch evidence.** Stage-6 residuals come from a smoothed track and are
strongly correlated in time, so they make every offset look "persistent".
Strength and persistence are therefore judged on the **epoch-by-epoch
Stage-4 position − TLE**, `stage4_d2_tle` (χ²₃ 0.9973 = 14.16):

- an epoch is *strong* when that value exceeds the limit;
- it is *persistent* when it is strong in ≥ `--min-persistent-epochs` (5)
  consecutive anomalous epochs.

**Station test** (epochs with ≥ 3 stations) uses two methods:

| Method | Rule |
|---|---|
| pairwise | pairs with `pair_separation_z > 4` are inconsistent. A station is picked if it belongs to every inconsistent pair and every pair without it is consistent. One inconsistent pair is `AMBIGUOUS`. |
| reference | used only with a precision TLE. The Stage-4 solution disagrees with the TLE, and removing exactly one station makes the disagreement go away. |

If both methods name the same station, that is cross-method agreement. If
only one names a station, that station is used. If they conflict, the
result is `AMBIGUOUS`.

`STATION_SPECIFIC` (run level) needs a station picked at ≥ 3 epochs **and**
at ≥ half the epochs with ≥ 3 stations.

**Tags** are non-exclusive:

| Tag | Source |
|---|---|
| `STATISTICAL_ANOMALY` | the anomaly statistic is above its threshold |
| `STATION_SPECIFIC` | the station test above |
| `CALIBRATION_CORRELATED` / `LOCALISATION_CORRELATED` / `TIMING_CORRELATED` | a supported correlation for a contributing station or the network. At the epoch, also `TIME_MISMATCH` (timing), or Stage-4 `HIGH_LINE_RESIDUAL`/`HIGH_CHI2` or an inconsistent pair (localisation) |
| `GEOMETRY_CORRELATED` | a supported geometry correlation |
| `TRAJECTORY_MODEL_CORRELATED` | Stage 5 raised `DYNAMICS_MISMATCH`/`HIGH_NIS` or restarted, or a NIS correlation is supported. At the epoch, `INNOVATION_OUTLIER`/`FILTER_DIVERGENCE` |
| `REFERENCE_LIMITED` | Stage 6 is not a precision validation, or the epoch is `TLE_STALE`/`TLE_REFERENCE_DISPERSION_HIGH` |
| `NETWORK_CONSISTENT` | Level-4 conditions below |
| `UNEXPLAINED` | Level 5 |
| `NORMAL` | no anomaly |

**Evidence level** (per epoch; overall = maximum). These are evidence
levels, not cause labels.

| Level | Rule |
|---|---|
| 0 | no statistical anomaly |
| 1 | statistical anomaly |
| 2 | + the station is picked by the station test, and the station is `STATION_SPECIFIC` |
| 3 | + calibration / timing / localisation evidence |
| 4 | `NETWORK_CONSISTENT`: ≥ 3 stations; every LOO shift within χ²_dof(0.9973), so no subset removes the anomaly; no inconsistent pair; no station picked; no instrument, geometry or trajectory-model tag; reference not limited |
| 5 | `UNEXPLAINED`: Level 4 **and** strong **and** persistent **and** calibration metadata available |

Two caps apply. `REFERENCE_LIMITED` caps the level at 3. Missing calibration
metadata caps it at 4, because missing evidence must not count as clean.

An epoch with fewer than 3 stations cannot have the subset test
(`SUBSET_TEST_UNAVAILABLE`), so it cannot reach Level 4.

**Confidence** (`LOW` / `MODERATE` / `HIGH`; never a probability) counts six
boolean components:

- sample support;
- effect size;
- persistence;
- cross-method agreement;
- reference quality;
- data completeness.

`HIGH` needs all six. `MODERATE` needs ≥ 4, with reference quality and data
completeness not both false. Anything less is `LOW`.

## 7. Outputs

```text
results/stage7/<YYYYMMDD>/NORAD_<id>/
├── diagnosis.{parquet,csv}          one row per Stage-6 epoch: residuals, D², Stage-4 − TLE D², strength/persistence,
│                                    NIS and flags, geometry (measured and predicted), largest line/direction residual
│                                    and LOO shift, station-test pick, tags, evidence level, confidence, quality flags
├── station_influence.{parquet,csv}  long: per (epoch, station) LOO solution, shift (3-D, RTN), loo_mahalanobis/dof,
│                                    degraded geometry, direction residuals, line residual, metadata, flags
├── pairwise.{parquet,csv}           per (epoch, station pair) diagnostic two-line solution
├── correlations.csv                 every test: scope, station, family, x, y, coefficients, p-values, n, n_eff, tag
├── summary.json                     overall level/tags/confidence, anomaly and unexplained intervals, station test,
│                                    station summaries, strongest correlation per family, reference and Stage-5 status
├── stage7_manifest.json             inputs and hashes, provenance chain, alignment, LOO/pairwise, correlation settings,
│                                    evidence rules with exact thresholds, policy, software, git, command line, outputs
├── diagnostics/  rtn_evidence  station_influence  station_direction_residuals  calibration_correlations
│                 timing_vs_along_track  geometry_vs_residual  influence_heatmap   (.png)
└── logs/stage7.log
```

**Station summaries** (§24) report, per station, over the diagnosed epochs:

- usable fraction;
- median direction and line residual;
- median LOO shift and LOO D²;
- zero-gain fraction;
- median flagged fraction;
- gain-solution changes;
- median SNR;
- visibility–image offset;
- median and MAD of the time offset;
- evidence tags and health flags.

No combined health score is produced. The optional score of §25 is not
implemented, so no opaque number exists.

Every table carries `norad_id` (and `station_id` where relevant). Runs for
several satellites can be concatenated without schema changes (§39).

**Exit status:**

| Code | Meaning |
|---|---|
| 0 | diagnosis completed. A limited result, such as `REFERENCE_LIMITED`, is still a completed diagnosis |
| 3 | fewer than 3 diagnosable epochs (insufficient data, e.g. an invalid upstream chain). No diagnosis is manufactured |
| 2 | invalid input, configuration or provenance chain |
| 1 | unexpected failure |

**Quality policy.**

- **`strict`** (default) diagnoses only Stage-6 `valid` rows that link to
  Stages 5 and 4.
- **`permissive`** also diagnoses rows with a finite residual, and tolerates
  hash mismatches (recorded).

## 8. Synthetic validation

`stage7.synthetic.build_chain` builds each scenario through the **real** chain
with genuine provenance and hashes:

1. Truth: the PRN-26 TLE trajectory by the independent TEME route of Stage 6,
   or a known modification of it.
2. Synthetic Stage-3 directions from the example stations (0.005° noise),
   with metadata columns.
3. A matching synthetic Stage-2 manifest with per-file calibration.
4. `tart-triangulate` → `tart-track` → `tart-validate` against a fresh TLE
   (0.35 d) → `tart-diagnose`.

`tests/stage7/` has 63 tests. The demo is `scripts/stage7_synthetic_demo.py`;
outputs are in [`examples/stage7`](examples/stage7), with the report in
[`demo_report.txt`](examples/stage7/demo_report.txt). **9 of 9 scenarios
match their expectation.**

| Scenario (spec) | Injected | Recovered tags | Level (confidence) | Station test / LOO |
|---|---|---|---|---|
| clean (35.1) | nothing | `NORMAL` | 0 (HIGH) | no picks |
| biased station (35.2, 36) | 0.036° bias at Madagascar | `STATISTICAL_ANOMALY`, `STATION_SPECIFIC`, `LOCALISATION_CORRELATED` | 3 (MODERATE) | Madagascar picked at 36/36 epochs; most influential: Madagascar. Removing Madagascar cuts the error vs truth below 0.3×; removing a healthy station does not |
| calibration (35.3) | the same bias only where the flagged fraction is high (40 % of epochs) | … `CALIBRATION_CORRELATED` (Madagascar `flagged_fraction` vs line residual, r = +0.83) … | 3 (MODERATE) | Madagascar picked at 12 epochs (below the run threshold) |
| timing (35.4) | Mauritius samples 1–5 s late (recorded), matched without interpolation | … `TIMING_CORRELATED` (largest time offset vs Stage-4 − TLE, r = +0.78) … | 3 (MODERATE) | Mauritius picked at 11 epochs |
| poor geometry (35.5) | blocks of 4-station and narrow 2-station (2.8°) epochs, no bias | `NORMAL`, `GEOMETRY_CORRELATED` (r = +0.72) | 0 (HIGH) | not network-consistent |
| dynamics (35.6) | unmodelled 0.5 m/s² acceleration | `STATISTICAL_ANOMALY`, `TRAJECTORY_MODEL_CORRELATED` (Stage 5: `HIGH_NIS`, `DYNAMICS_MISMATCH`) | 1 (HIGH) | — |
| stale TLE (35.7) | 10 s offset against a 20-day-old TLE | `STATISTICAL_ANOMALY`, `REFERENCE_LIMITED` | 1, capped (MODERATE) | — |
| network offset (35.8/35.9) | 10 s network-wide time shift (39 km along-track) | `STATISTICAL_ANOMALY`, `NETWORK_CONSISTENT`, `UNEXPLAINED` | **5 (HIGH)** | no picks; LOO D² median < 5 |
| missing Stage 2 (§33) | same, no calibration metadata | `STATISTICAL_ANOMALY`, `NETWORK_CONSISTENT` + `CALIBRATION_METADATA_MISSING` | 4, Level 5 blocked (MODERATE) | — |

Further checks in the tests:

- With `--min-persistent-epochs 100`, the network case stops at Level 4.
- No physical-cause label appears in any output.
- Upstream files are unchanged.

**Timing signatures.** A network-common time offset appears as an
along-track residual: in the network case |dT| > 5|dR| and > 5|dN| at more
than 90 % of epochs. A *single* station's offset is different. It maps
through the geometry mostly along that station's line of sight: in the
timing scenario the Stage-4 − TLE error correlates with the offset in all
three components, but its RMS is 39 km radial against 3 km along-track.
"Timing appears along-track" therefore holds only for a network-wide offset,
and Stage 7 does not assume it.

**Power.** The geometry case is detected in 5 of 6 noise realisations
(partial ρ 0.41–0.72); the miss falls just below the effect-size threshold.
Magnitudes of residuals drawn at two different scales overlap heavily, so
Spearman correlations between geometry and |residual| stay moderate even for
a large σ contrast.

## 9. Missing data

Missing diagnostics never crash the analysis, and missing evidence is never
treated as clean:

- **No Stage-2 manifest:** `CALIBRATION_METADATA_MISSING`, confidence
  lowered, Level 5 blocked.
- **No Stage-3 row** for a contributing station: `STAGE3_LINK_MISSING`, and
  that station is not re-triangulated.
- **Fewer than 3 stations:** `SUBSET_TEST_UNAVAILABLE`.

## 10. Real-data run (2026-06-25, GPS PRN 26, NORAD 40534)

Outputs are in [`examples/stage7_real`](examples/stage7_real).

**Inputs.** The whole chain verifies: the Stage-6, Stage-5, Stage-4 and
Stage-3 example manifests, and the Stage-2 manifest in `../sat_position`,
whose hash matches. The full re-triangulation reproduces Stage 4 to 0 m.

**Spec form on the two-station run (`diagnostic_min2`): exit 0.**

| | |
|---|---|
| diagnosed epochs | 36 of 69 (the others have no Stage-5 state) |
| anomalous | 36 (threshold 21.7, conditional statistic) |
| evidence level | **3, `REFERENCE_LIMITED` (MODERATE)**: 26 epochs at Level 3, 10 at Level 1 |
| tags | `STATISTICAL_ANOMALY`, `LOCALISATION_CORRELATED` (11 of 37 station pairs have `pair_separation_z > 4`; Stage-4 `HIGH_CHI2`), `TRAJECTORY_MODEL_CORRELATED` (Stage 5 `DYNAMICS_MISMATCH`), `REFERENCE_LIMITED` (103-day back-propagated TLE) |
| leave-one-out | only **1** epoch has three stations, so the most influential station is **undetermined** and the station test cannot run |
| correlations | 130 tests, **none supported**: 59 lack samples, 30 lack variation, 39 are not supported, 2 are expected. The strongest localisation candidate is Mauritius visibility–image separation vs along-motion residual (raw 0.83, partial 0.58), but with n_eff ≈ 17 its Holm p is 0.27, so it is not supported |

Station summaries (diagnosed epochs):

| Station | Epochs | Median direction residual | Median line residual | Median SNR | Flagged | Zero-gain fraction |
|---|---|---|---|---|---|---|
| South Africa | 36 | 0.27° | 117 km | 10.0 | 0.62 | 0.375 (9 of 24 antennas) |
| Mauritius | 35 | 0.31° | 127 km | 9.5 | 0 | 0 |
| Madagascar | 1 | 1.42° | 432 km | 9.7 | 0 | 0 |

All three carry the Stage-3 `CALIBRATION_WARNING`. Madagascar's single
contributing direction is 1.4° from the network trajectory, consistent with
the other source Stage 3 found there. Mauritius and South Africa wander
slowly by 0.1–0.6° (see `station_direction_residuals.png`).

**Default Stage-6 run: exit 3.** There are no diagnosable epochs (no valid
Stage-5 rows), and no diagnosis is manufactured.

**Interpretation.** The real result is *limited*, not anomalous in any strong
sense:

- the reference is a stale TLE;
- the trajectory model is inconsistent with the data;
- localisation is internally inconsistent;
- a three-station test is possible at a single epoch.

None of this supports a station-specific or network-level attribution.

Two things would strengthen it:

- a contemporaneous TLE, so that `REFERENCE_LIMITED` lifts and the reference
  method of the station test becomes available;
- more three-station epochs, i.e. a cleaner MG-Tana track.

## 11. Limitations

- The evidence rules are deliberately simple and transparent. Thresholds are
  in the manifest; they are not learned.
- Correlations are within one pass; the effective-sample-size and
  time-detrending corrections are approximations. Smooth, non-linear common
  trends can still produce rare false positives.
- With three stations the subset test is weak (§3), and with two it is
  unavailable.
- Direction residuals use the Stage-5 trajectory, which the same measurements
  produced. A strong single-station bias pulls that trajectory and is
  partially absorbed.

## 12. Command name

Stage 7 adds `tart-diagnose` (`tart_gnss.stage7.cli:main`). The public
commands are:

- `tart-stations-update`
- `tart-stations-list`
- `tart-common-gnss`
- `tart-retrieve`
- `tart-localise`
- `tart-triangulate`
- `tart-track`
- `tart-validate`
- `tart-diagnose`

There are no aliases. `tests/test_cli_names.py` checks that the entry points
match this set exactly, and that every `prog=` matches its entry point.
