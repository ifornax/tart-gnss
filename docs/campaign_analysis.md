# Campaign analysis across completed runs (`tart-campaign`)

The single-pass pipeline (Stages 1–7) asks whether the network can
reconstruct and diagnose *one* pass. The campaign layer asks how well it
performs **statistically**, across satellites, passes, stations, geometries,
constellations and observing conditions.

It is not another processing stage. `tart-campaign` reads **completed**
Stage-7 runs, plus the upstream products they recorded. It never runs
`tart-localise`, `tart-triangulate`, `tart-track`, `tart-validate` or
`tart-diagnose`; a test makes every stage entry point raise during a
campaign. An incomplete run is reported as incomplete, never finished
silently.

Specification: [`tart_gnss_campaign_analysis_implementation.md`](../tart_gnss_campaign_analysis_implementation.md).

## 1. Inputs, provenance and the campaign unit

**Unit.** One completed Stage-7 run = one satellite × one common observing
window. Each run gets a `run_id` of the form `<YYYYMMDD>_NORAD<id>_<HHMMSS>`,
from its first diagnosed epoch. A numeric suffix makes the ID unique when two
analyses cover the same window. The original manifest path and SHA-256 are
kept.

**Discovery.** Runs come from `--input-root`, a recursive search for files
named exactly `stage7_manifest.json` (repeatable), and/or from explicit
`--stage7-manifest` arguments. Runs are never inferred from CSV or Parquet
files. A manifest found twice is listed once.

**Loading.** Each run is checked as follows:

- `stage == 7`, with a schema version in the supported set (currently 1). A
  future version fails as `SCHEMA_MISMATCH`.
- Its own outputs must be present: `summary.json`, `diagnosis` and
  `station_influence`.
- Upstream products are reached only through the `provenance_chain` that
  Stage 7 recorded. Paths resolve as recorded, then relative to the Stage-7
  run's working directory, and every upstream manifest is **re-hashed**
  against the hash Stage 7 recorded.

The products read are:

| Stage | Products |
|---|---|
| 6 | residuals and summary |
| 5 | trajectory and manifest |
| 4 | positions and station contributions |
| 3 | manifest (satellite metadata) |

A missing, changed or mismatched link makes the run `BROKEN_PROVENANCE`.

Per run the campaign records:

- path, SHA-256, git commit and pipeline version;
- NORAD ID, satellite name and constellation;
- date, start, end and duration;
- station set;
- TLE source and epoch, and TLE age;
- precision validation;
- overall evidence level, tags and confidence;
- the most influential station.

**Constellation** comes from structured metadata only: the Stage-3 manifest
`satellite.constellation`, which is copied from Stages 1–2. The canonical
values are GPS, Galileo, BeiDou, GLONASS and Unknown. Names are never
parsed. A missing value is `Unknown`, flagged `MISSING_CONSTELLATION_METADATA`.

## 2. Policy, filters, cohorts

Every run is either included or listed in `excluded_runs.csv` with a reason.
Runs are never dropped silently.

| Reason | When |
|---|---|
| `SCHEMA_MISMATCH` | unsupported Stage-7 schema |
| `MISSING_STAGE7_OUTPUT` | unreadable manifest, or missing Stage-7 tables |
| `BROKEN_PROVENANCE` | `strict` only (the default). Under `permissive` the run is kept and flagged |
| `INSUFFICIENT_EPOCHS` | fewer than `--min-epochs-per-run` (3) diagnosed epochs |
| `NO_PRECISION_VALIDATION` | with `--precision-validation-only` |
| `FILTERED_*` | `--norad-id`, `--constellation`, `--station`, `--start`/`--end`, `--min-elevation`, `--max-tle-age-days`, `--max-evidence-level` |
| `DUPLICATE_RUN` | identical Stage-7 manifest content already included |

Included runs fall into one of two **cohorts**:

- `precision_validation`: the Stage-6 precision validation is true;
- `descriptive_only`.

Scientific summaries of TLE-relative residuals use the precision cohort by
default. `--no-precision-only` analyses all included runs, and
`--precision-validation-only` excludes non-precision runs altogether. The
counts of both cohorts are always reported.

## 3. Terminology

The campaign never calls an observed − TLE residual an absolute accuracy:

| Term | Meaning |
|---|---|
| **measurement uncertainty** | formal 1σ of the Stage-3/4/5 fits |
| **TLE-relative residual** | Stage-6 observed − TLE. A performance measure only in the precision cohort |
| **trajectory residual** | Stage-4 positions about the Stage-5 trajectory (RMS); not a TLE comparison |
| **network repeatability** | pass-to-pass scatter of per-pass summaries |

## 4. Statistics: runs are the independent units

Epochs within a pass are strongly correlated: Stage-6 residual lag-1
autocorrelation is about 0.9. Every inference therefore aggregates per run
first:

- population statistics use one per-run median per run;
- grouped comparisons use per-run values;
- dependence analyses use per-run values;
- **confidence intervals come from a run-level (cluster) bootstrap**:
  resample runs with replacement, percentile 95 %, `--bootstrap` (1000)
  resamples, fixed `--seed`.

Pooled-epoch distributions are written for display only, and never yield a
p-value. A test shows that resampling by run gives intervals more than 5×
wider than a naive epoch bootstrap on correlated epochs. Bootstrap coverage
of the median of a known lognormal population is 0.88–0.99 over 200
replications (tested).

**Comparisons** are made only when every group has ≥ `--min-runs-per-group`
(3) runs; smaller groups are flagged `under_sampled`.

| Groups | Test | Effect size |
|---|---|---|
| 2 | Mann-Whitney U | rank-biserial correlation, with its sign stated, plus the difference in medians |
| > 2 | Kruskal-Wallis | ε² |

**Dependence** (geometry, elevation, SNR, TLE age) uses the Spearman
correlation and the Theil-Sen slope with its 95 % interval, on per-run
values. No machine learning is involved.

**P-values are permutation p-values.** With a handful of runs and tied
per-run values, large-sample approximations can be orders of magnitude too
small. For 3 against 13 tied runs, asymptotic Mann-Whitney gives about 10⁻⁴;
the exact value is 1/560. Group-label permutations, or pairings for Spearman,
are therefore enumerated exactly when there are at most 10⁵ of them (Spearman:
n ≤ 8). Otherwise 10⁴ seeded random permutations are used, and `p_method`
records which. With fewer than 5 runs, a dependence is flagged because its
Theil-Sen interval is not meaningful.

**Fractions** of runs carry Wilson 95 % intervals. Epoch-level and
station-level fractions get no interval, because their rows are correlated
within a run.

**Sample-size guidance:**

| Runs in the analysis cohort | Tier |
|---|---|
| < 3 | `descriptive_only`, exit 3, `INSUFFICIENT_CAMPAIGN_SIZE` |
| 3–9 | `preliminary` |
| ≥ 10 | `basic_inference` |
| ≥ 30 | `grouped_comparisons` |

## 5. What is computed

| Topic | Content |
|---|---|
| direction (Stage 3) | median σ_az, σ_el and angular σ; visibility–image disagreement; direction residual vs the Stage-5 trajectory, kept separate from the formal σ |
| 3-D localisation (Stage 4) | 3-D σ and the **principal axes** of the position covariance (`median_sigma_major/intermediate/minor_m`), since the uncertainty is strongly anisotropic |
| trajectory (Stage 5) | filtered and smoothed σ, velocity σ, trajectory residual vs Stage 4, median NIS, gate-exceedance fraction, lag-1 innovation autocorrelation, filter resets, dynamics-mismatch fraction |
| TLE-relative residual (Stage 6) | 3-D, radial, along-track, cross-track and velocity residual: median, mean, RMS, MAD, 16/84/95th percentiles, with bootstrap CIs (`--reference-frame rtn` or `gcrs`) |
| TLE age | per-run residuals in \|age\| bins (0–1, 1–3, 3–7, 7–14, > 14 d) with n, median, RMS and MAD; Spearman/Theil-Sen against \|age\| and signed age. Stale runs appear as descriptive; no causation is inferred |
| elevation, SNR | run × station metrics in 10° elevation bins and SNR quantile bins (linear: SNR spans less than a decade) |
| geometry | per-run metrics in crossing-angle bins; σ vs crossing angle **within each station count**, since station count dominates; σ by number of stations; leave-one-out shift vs crossing angle |
| station stability | per station across runs: runs, epochs, usable fraction, median direction residual and its run-to-run **standard deviation**, MAD and 84th percentile, line residual, LOO influence, flagged fraction, zero-gain and gain-change incidence, timing median and MAD, and the fraction of runs tagged `STATION_SPECIFIC` / `CALIBRATION_` / `TIMING_` / `LOCALISATION_CORRELATED` |
| comparisons | constellation, satellite (NORAD) and station-set summaries; Mann-Whitney / Kruskal-Wallis on SNR, direction σ, Stage-4 σ and residual |
| evidence | level fractions 0–5 by **run, epoch and station separately**; tag frequencies; tag rates against station, constellation, crossing-angle bin and TLE-age bin |
| repeatability | per satellite with ≥ 2 passes: between-pass std and MAD of the per-pass median R/T/N against the within-pass median σ (with their ratio), and the between-pass scatter of each station's direction bias and influence. Evidence-tag repeatability is the mean Jaccard similarity of tag sets |
| balance | share of the largest satellite, constellation, station set, day and elevation range |

Stability uses the run-to-run standard deviation rather than the MAD. A robust
MAD hides a station that is biased in only some runs, which is exactly the
instability of interest; a test checks this.

No combined "best station" score is produced; the station table is a
transparent comparison. No physical-cause rates are produced: no
manoeuvre, clock-failure, hardware-failure or RFI rates. The reported
quantities are tag fractions such as `NETWORK_CONSISTENT`, `UNEXPLAINED`,
`STATION_SPECIFIC` and `REFERENCE_LIMITED`.

**Campaign quality flags:**

| Flag | Condition |
|---|---|
| `INSUFFICIENT_CAMPAIGN_SIZE` | < 3 runs in the analysis cohort |
| `INSUFFICIENT_GROUP_SIZE` | a constellation group is under-sampled |
| `LOW_PRECISION_VALIDATION_FRACTION` | < 0.5 |
| `BROKEN_PROVENANCE_PRESENT` | |
| `DOMINATED_BY_ONE_SATELLITE` | > 50 % of runs |
| `DOMINATED_BY_ONE_STATION_SET` | > 80 % of runs |
| `TLE_AGE_IMBALANCE` | > 80 % of runs in one age bin |
| `MISSING_CONSTELLATION_METADATA` | |

Otherwise the flag is `GOOD`.

## 6. Outputs

```text
results/campaign/<campaign_name>/
├── campaign_manifest.json        inputs + hashes, included/excluded runs, filters, grouping, bootstrap (n, seed),
│                                 policy, software, git commit, command line, outputs
├── campaign_runs.{parquet,csv}   one row per run (section 1 fields and all per-run performance metrics)
├── campaign_epochs.{parquet,csv} one row per diagnosed epoch
├── campaign_stations.{parquet,csv}  one row per run x station
├── excluded_runs.csv             run_id, manifest_path, norad_id, reason, quality_flags, detail
├── summary.json                  counts, span, cohort sizes, bootstrap-CI medians, performance, evidence and tag
│                                 fractions, station highlights, comparisons, dependence, repeatability, balance, flags
├── tables/  station_summary  constellation_summary  satellite_summary  evidence_summary  tle_age_summary
│            geometry_summary  elevation_summary  snr_summary  station_set_summary  anomaly_rates  grouped_by_<field>
├── diagnostics/  residual_distribution  rtn_distribution  residual_vs_elevation  residual_vs_snr
│                 uncertainty_vs_crossing_angle  station_stability  constellation_comparison
│                 evidence_level_distribution  tle_age_dependence  repeatability   (.png; + .pdf with --publication-plots)
└── logs/campaign.log
```

`--publication-plots` gives larger fonts, 200-dpi PNG plus vector PDF, and
labelled axes with units.

**Exit status:**

| Code | Meaning |
|---|---|
| 0 | campaign analysis completed (≥ 3 runs in the analysis cohort) |
| 3 | completed, but too few runs for campaign inference; descriptive outputs are still written |
| 2 | invalid input or configuration |
| 1 | unexpected failure |

## 7. Synthetic campaign validation

`tart_gnss.campaign.synthetic` builds a 17-run campaign through the **real**
Stage 3–7 chain, with known population properties:

- **Satellites:** three. The PRN-26 example (GPS, 0.005°, SNR 20), plus two
  relabelled neighbours 2° away in mean anomaly: SYN-GAL (Galileo, NORAD
  90001, 0.008°, SNR 12) and SYN-BDS (BeiDou, 90002, 0.006°, SNR 16).
- **Repeated passes:** on later sidereal days, when the ground track repeats.
- **Unstable station:** Madagascar is biased by 0.036° in three runs.
- **Geometry:** two-station runs on each station pair (crossing angles about
  2.8, 6.6 and 8.6°) and a four-station run including Namibia.
- **TLE age:** each GPS run carries an along-track drift of 0.4 s per day of
  TLE age (about 1.55 km/day). Ages run from 0.35 to 20 d, and the 20-day run
  is stale, i.e. descriptive only.
- **Repeatability:** three SYN-GAL passes offset by −1, 0 and +1 s along
  track, a known between-pass standard deviation of |v| × 1 s = 3.87 km.

`tests/campaign/` has 58 tests: one per module, plus §47.1–47.6 on this
campaign. The demonstration is `scripts/campaign_synthetic_demo.py`; outputs
are in [`examples/campaign`](examples/campaign), with the report in
[`demo_report.txt`](examples/campaign/demo_report.txt). It has 17 runs, 3
satellites, 3 constellations and 4 stations, and its tier is
`basic_inference`. 16 runs are precision validations (fraction 0.94) and 1
is descriptive.

**Measurement uncertainty** (per-run medians, bootstrap 95 % CI):

| Quantity | Value |
|---|---|
| Stage-4 3-D σ | 16.6 km (16.6–20.2) |
| Stage-4 major / minor axis | 16.5 / 1.09 km, about 15:1 anisotropy along the line of sight |
| Stage-5 smoothed 3-D σ | 4.47 km (4.46–5.43) |

**TLE-relative residual** (precision cohort):

| Component | Value |
|---|---|
| 3-D | 13.8 km (8.1–18.4) |
| radial | 5.8 km |
| along-track | 3.3 km |
| cross-track | −0.4 km |

These include the injected age drift; nothing here is an absolute accuracy.

**Recovery of the injected properties** (§47):

| Property | Recovered |
|---|---|
| 47.1 clean | clean runs at Level 0; a fresh-TLE GPS pass has a residual within a few Stage-5 σ |
| 47.2 unstable station | Madagascar is the only station tagged `STATION_SPECIFIC` (18 % of runs), the most frequently flagged, with the largest run-to-run std (0.0076°) and 84th percentile (0.019°) of its direction residual. Its LOO influence is *not* higher: with three stations, removing any station moves a biased solution (see the Stage-7 doc) |
| 47.3 geometry | Stage-4 σ is 22.9 / 16.6 / 12.2 km for 2 / 3 / 4 stations. Among two-station runs, σ falls with crossing angle (Spearman −1, Theil-Sen −0.085 dex/deg). With only 3 runs the exact p is 1/3, so this is suggestive only, and the dependence is flagged |
| 47.4 TLE age | Spearman +0.59 (permutation p 0.014); Theil-Sen slope +1.13 km/day (95 % 0.12–1.53), against ~1.55 km/day injected and diluted by the drift-free runs |
| 47.5 constellation | Galileo SNR is below GPS: Mann-Whitney exact permutation p = 0.0018 (= 1/560; asymptotic would give 0.00016), rank-biserial −1.00. BeiDou (1 run) is `under_sampled` |
| 47.6 repeatability | SYN-GAL between-pass along-track std is 3.58 km, against 3.87 km known (7 % low with 3 passes), and 7.5× its within-pass σ (0.48 km) |

The run-level evidence fractions are:

| Level | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| fraction of runs | 0.24 | 0.12 | 0.18 | 0.12 | 0.35 | 0 |

Level 5 is absent because the drifting GPS runs are network-consistent but
the evidence is capped. Tags: `NETWORK_CONSISTENT` 0.35 [0.17, 0.59],
`STATION_SPECIFIC` 0.18 [0.06, 0.41], `REFERENCE_LIMITED` 0.06.

**One finding from building this campaign.** When a network-wide offset and a
station bias occur together, Stage 7's TLE-based station test and its pairwise
test can name different stations. Stage 7 then reports `AMBIGUOUS` and no
`STATION_SPECIFIC` tag. That is a conservative, documented limitation, so the
synthetic unstable-station runs carry no age drift.

## 8. Real campaign (available runs)

```bash
tart-campaign --input-root docs/examples/stage7_real --campaign-name real_tart_gnss --output-dir docs/examples/campaign_real
```

Outputs are in [`examples/campaign_real`](examples/campaign_real). The run
**exits 3** (`INSUFFICIENT_CAMPAIGN_SIZE`, `LOW_PRECISION_VALIDATION_FRACTION`,
`INSUFFICIENT_GROUP_SIZE`):

| | |
|---|---|
| discovered | 2 Stage-7 runs (2026-06-25, PRN 26) |
| included | 1 (the two-station diagnostic run); the default run is excluded as `INSUFFICIENT_EPOCHS` (0 diagnosed epochs) |
| precision-valid | 0 (TLE 103 days old) |
| satellites / stations / constellations | 1 / 3 / GPS |
| span | 2026-06-25 09:18:35–09:24:22 (347 s) |
| descriptive | Stage-4 σ 547 km (major axis 546 km, minor 26 km); Stage-5 σ 130 km; evidence Level 3 |

With one non-precision run, there are no population results and none are
claimed. The outputs exist to show the tool runs on real data.

## 9. Command name

`tart-campaign` (`tart_gnss.campaign.cli:main`). The full workflow is:

```text
per pass:  tart-common-gnss → tart-retrieve → tart-localise → tart-triangulate → tart-track → tart-validate → tart-diagnose
campaign:  tart-campaign
```

`tests/test_cli_names.py` checks that the entry points match this set of
`tart-*` names exactly.
