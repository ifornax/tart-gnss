# tart-gnss

**A pipeline for localising GNSS satellites with the TART radio-telescope network.
Stage 1 finds satellites that three TART stations can see at the same time.
Stage 2 retrieves, validates and synchronises the real TART observations for
those windows. Stage 3 measures the satellite's direction at each station
from the complex visibilities. Stage 4 intersects those measured lines of
sight to give the satellite's 3-D position with its covariance. Stage 5
turns those positions into a continuous position/velocity trajectory with a
full covariance, without using any TLE. Stage 6 then compares that trajectory
with a TLE prediction, as observed − TLE residuals. Stage 7 traces the whole
chain back to Stage 2 and reports evidence about where a discrepancy comes
from, without asserting a physical cause. The campaign layer turns many
completed runs into population statistics.**

| Stage | Command | What it does |
|---|---|---|
| 1 | `tart-common-gnss` | Propagates GNSS orbits from TLEs (SGP4 via Skyfield) and finds every continuous window in which a satellite is above 20° at three stations at once. Windows are refined to 0.1 s, ranked, and stored with ITRS positions and TLE provenance. |
| 1→2 | `tart-archive-scan` | Utility between Stages 1 and 2. Lists the TART archive once per station and UTC day (with a checksummed cache), turns the filenames into conservative per-station *candidate* coverage, intersects the stations, matches the shared segments with **every** Stage-1 window, and ranks the in-band candidates. Each candidate gets a safe `tart-retrieve --dry-run` command. With `--all-stations` it scans the whole registry for intervals and triads at ≥ 3 distinct physical sites. Nothing is downloaded; results are archive-coverage candidates, not validated data. |
| 2 | `tart-retrieve` | Takes one Stage-1 window, queries the TART visibility archive, downloads the native HDF5 files with SHA-256 provenance, validates timestamps, frequency, station identity and calibration, measures the real three-station overlap, converts each station to its own Measurement Set with `tart2ms`, and writes `manifest.json` for Stage 3. |
| 3 | `tart-localise` | Reads the Stage-2 manifest and, for each station independently, fits a point-source model directly to the complex visibilities in time bins (profiled-χ² grid search + robust `least_squares`), giving Az/El with covariance, SNR, goodness-of-fit and quality flags, cross-checked by a dirty-image peak finder. Writes CSV/Parquet tracks, `stage3_manifest.json` and diagnostic plots for Stage 4. |
| 4 | `tart-triangulate` | Reads the Stage-3 directions and the frozen station coordinates, matches stations in time (nearest within tolerance, no interpolation), and intersects the ECEF sight lines — covariance-weighted (range-scaled Stage-3 angular errors) or closed-form unweighted — giving x, y, z (ITRS/ECEF) with a 3×3 covariance, station ranges, line residuals, geometry and timing diagnostics, quality flags and an optional bootstrap. TLEs never enter the fit. |
| 5 | `tart-track` | Reads the Stage-4 positions and full 3×3 covariances, rotates them ITRS/ECEF → GCRS (Skyfield), initialises a position/velocity state from the first Stage-4 measurements only, and runs an extended Kalman filter (two-body by default; constant-velocity and J2 available) with χ²₃ innovation gating, followed by an RTS smoother. Writes the filtered and smoothed GCRS state, the ECEF position, the full 6×6 covariances, innovations, NIS/whiteness diagnostics, quality flags, `stage5_manifest.json` and seven plots. TLEs are never used. |
| 6 | `tart-validate` | Propagates exactly one selected TLE (`--tle-file`, or `--tle-source stage1`/`download`) with SGP4 at the exact Stage-5 epochs and transforms it TEME → GCRS. Reports **observed − TLE** residuals: Cartesian, radial/along-track/cross-track (RTN), angular, range and an along-track time equivalent. Also gives the Stage-5 covariance rotated into RTN, measurement-conditional Mahalanobis distances (and reference-aware ones with an ensemble of ≥ 3 TLEs), TLE age, and a precision-validation status. Writes `residuals`, `tle_trajectory`, `summary.json`, `stage6_manifest.json` and eight plots. The TLE is never fitted, and the Stage-5 run is never modified. |
| 7 | `tart-diagnose` | Follows the recorded provenance from Stage 6 back to Stage 2, verifying the hashes, and aligns everything by time and station. Recomputes leave-one-station-out and pairwise solutions with the Stage-4 solver, compares every station direction with the Stage-5 trajectory, and tests calibration, localisation, timing, geometry and trajectory-model correlations (time-detrended, autocorrelation-aware, effect size plus Holm). Assigns non-exclusive evidence tags, rule-based evidence levels 0–5 and a descriptive confidence. Writes `diagnosis`, `station_influence`, `pairwise`, `correlations.csv`, `summary.json`, `stage7_manifest.json` and seven plots. No physical cause is asserted, and no upstream output is modified. |
| campaign | `tart-campaign` | Aggregates many completed Stage-7 runs, reading only their recorded upstream products (re-hashed; nothing is rerun). Produces run, epoch and station population tables, and reports: measurement uncertainty (including the Stage-4 covariance principal axes); TLE-relative residuals (precision cohort); station stability; constellation, satellite and station-set comparisons (with effect sizes and minimum group sizes); geometry, elevation, SNR and TLE-age dependence; evidence-level and tag rates; and pass-to-pass repeatability. Confidence intervals come from a run-level bootstrap, and runs — never epochs — are the independent units. |

> **Status:** Stages 1–7, the campaign layer and archive-aware candidate discovery complete, with dynamic TART station-registry integration (v0.10.0). A separate research track validates direction recovery on finite-range [tab-sim](https://github.com/chrisfinlay/tab-sim) simulations ([below](#simulation-validation-tab-sim-research)). The repository is private and under active development.

---

## Contents

- [Project context](#project-context)
- [Installation](#installation)
- [Stations](#stations)
- [Stage 1: common visibility (`tart-common-gnss`)](#stage-1-common-visibility-tart-common-gnss)
- [Archive scan: finding archived overlaps (`tart-archive-scan`)](#archive-scan-finding-archived-overlaps-tart-archive-scan)
- [Stage 2: observation retrieval (`tart-retrieve`)](#stage-2-observation-retrieval-tart-retrieve)
- [Stage 3: direction estimation (`tart-localise`)](#stage-3-direction-estimation-tart-localise)
- [Stage 4: 3-D localisation (`tart-triangulate`)](#stage-4-3-d-localisation-tart-triangulate)
- [Stage 5: trajectory tracking (`tart-track`)](#stage-5-trajectory-tracking-tart-track)
- [Stage 6: TLE validation (`tart-validate`)](#stage-6-tle-validation-tart-validate)
- [Stage 7: network diagnosis (`tart-diagnose`)](#stage-7-network-diagnosis-tart-diagnose)
- [Campaign analysis (`tart-campaign`)](#campaign-analysis-tart-campaign)
- [Simulation validation (tab-sim, research)](#simulation-validation-tab-sim-research)
- [TLEs and reproducibility](#tles-and-reproducibility)
- [Example results](#example-results)
- [Caveats](#caveats)
- [Repository layout](#repository-layout)
- [Development and testing](#development-and-testing)
- [Technical report](#technical-report)
- [Authors and licence](#authors-and-licence)

---

## Project context

The [Transient Array Radio Telescope (TART)](https://tart.elec.ac.nz) is a
low-cost, open-source 24-element all-sky interferometer operating at the GPS L1
frequency (1575.42 MHz, 2.5 MHz bandwidth). Three TART stations that image the
same satellite at the same time give three lines of sight, so the satellite's
3-D position can be inferred and compared with its orbital prediction.

| Stage | Task | Status |
|---|---|---|
| **1** | **Find satellites simultaneously visible from three TART stations** | **✅ `tart-common-gnss`** |
| **2** | **Retrieve, validate and synchronise the TART data; one Measurement Set per station** | **✅ `tart-retrieve`** |
| **3** | **Measure the satellite direction at each station from the visibilities (with image-domain cross-check)** | **✅ `tart-localise` — implemented** |
| **4** | **Multi-station 3-D localisation from the measured lines of sight** | **✅ `tart-triangulate` — implemented** |
| **5** | **TLE-independent Bayesian trajectory (position, velocity, 6×6 covariance) from the Stage-4 positions** | **✅ `tart-track` — implemented** |
| **6** | **Compare the measured trajectory with the TLE prediction (observed − TLE, Cartesian and RTN, with uncertainty context)** | **✅ `tart-validate` — implemented** |
| **7** | **Diagnose where a discrepancy comes from (station, instrument, timing, geometry, trajectory model, reference) as evidence levels** | **✅ `tart-diagnose` — implemented** |
| **campaign** | **Population statistics across many completed runs** | **✅ `tart-campaign` — implemented** |

Three conditions are kept strictly separate: **geometric visibility**
(Stage 1), **data availability** (Stage 2) and **radio detection** (Stage 3).
None is inferred from another.

Specifications: [`sat_localisation.md`](sat_localisation.md) (Stage 1),
[`tart_gnss_stage2_implementation.md`](tart_gnss_stage2_implementation.md) (Stage 2) and
[`tart_gnss_stage3_implementation.md`](tart_gnss_stage3_implementation.md) (Stage 3; method and
conventions in [`docs/stage3_direction_estimation.md`](docs/stage3_direction_estimation.md)) and
[`tart_gnss_stage4_implementation.md`](tart_gnss_stage4_implementation.md) (Stage 4; method in
[`docs/stage4_triangulation.md`](docs/stage4_triangulation.md)) and
[`tart_gnss_stage5_implementation.md`](tart_gnss_stage5_implementation.md) (Stage 5; method in
[`docs/stage5_tracking.md`](docs/stage5_tracking.md)) and
[`tart_gnss_stage6_implementation.md`](tart_gnss_stage6_implementation.md) (Stage 6; method in
[`docs/stage6_validation.md`](docs/stage6_validation.md)) and
[`tart_gnss_stage7_implementation.md`](tart_gnss_stage7_implementation.md) (Stage 7; method in
[`docs/stage7_diagnosis.md`](docs/stage7_diagnosis.md)) and
[`tart_gnss_campaign_analysis_implementation.md`](tart_gnss_campaign_analysis_implementation.md) (campaign layer;
method in [`docs/campaign_analysis.md`](docs/campaign_analysis.md)).

**Command names.** Every public command starts with `tart-` (`tart-stations-update`,
`tart-stations-list`, `tart-common-gnss`, `tart-retrieve`, `tart-localise`,
`tart-triangulate`, `tart-track`, `tart-validate`, `tart-diagnose`, `tart-campaign`). This is a stable public API, and future stages must follow it;
`tests/test_cli_names.py` enforces it.

## Installation

Requires Python ≥ 3.9 for Stage 1. Stage-2 Measurement-Set conversion needs
`tart2ms`, which requires Python 3.11–3.13.

```bash
git clone git@github.com:ifornax/tart-gnss.git
cd tart-gnss
pip install .                  # Stage 1 only (numpy, skyfield, sgp4, PyYAML)
pip install ".[stage2]"        # + archive access and HDF5 reading (minio, h5py)
pip install ".[stage2,ms]"     # + Measurement-Set conversion (tart2ms, python-casacore)
pip install -e ".[dev,ms]"     # development: everything plus pytest and matplotlib
```

Because the repository is private, a direct install needs SSH access:

```bash
pip install "tart-gnss[stage2,ms] @ git+ssh://git@github.com/ifornax/tart-gnss.git"
```

**After upgrading Python to 3.12 or later in an existing environment.** If
`pip` fails with `AttributeError: module 'pkgutil' has no attribute
'ImpImporter'`, the environment still has an old `pip` (before 23.1) or
`setuptools`, written for Python ≤ 3.11. `pip` cannot upgrade itself while
it crashes, so bootstrap it first, then reinstall this package (this also
rewrites the `tart-*` command scripts that still point at the old Python):

```bash
python -m ensurepip --upgrade
python -m pip install -U pip setuptools wheel
python -m pip install -e ".[dev,ms]"
```

## Stations

Coordinates are never hard-coded. Station metadata comes from the official
TART network registry, `https://map.elec.ac.nz/api/v1/telescopes` (its JSON
API; the map HTML is never scraped), through three files in
[`src/tart_gnss/config/`](src/tart_gnss/config):

| File | Maintained | Contents |
|---|---|---|
| `stations_registry.yaml` | **generated** by `tart-stations-update` | Default catalogue. One entry per valid registry telescope, keyed by canonical TART id, plus a `metadata` block (registry URL, retrieval time, SHA-256 of the raw snapshot, skipped entries, changes from the previous catalogue). |
| `station_aliases.yaml` | by hand | Convenience aliases: `mauritius → mu-udm`, `madagascar → mg-tana`, `south_africa → za-rhodes`, … |
| `stations.yaml` | by hand (legacy) | The original alias-keyed catalogue. Still accepted via `--station-file`. |

The raw registry payload behind the packaged catalogue is archived in
[`data/stations/`](data/stations) as `registry_<UTC>.json` with a
`.meta.json` sidecar (URL, retrieval time, SHA-256, station count) and a
`.diff.json` change record. Snapshots are never overwritten.

Identifiers resolve in this order: exact station id, normalised id
(case-insensitive, `-` = `_`), alias, then `tart_id`. So `za-rhodes`,
`za_rhodes` and `south_africa` select the same station.

```bash
tart-stations-list                      # table of the packaged catalogue, with aliases
tart-stations-list --archive-only --format json
tart-stations-update                    # fetch → data/stations/ snapshot + stations_registry.yaml
tart-stations-update --verify-archives  # also check <tart_id>/vis/ in the archive bucket
tart-stations-update --diff --dry-run   # show what changed, write nothing
tart-stations-update --verify-archives --output src/tart_gnss/config/stations_registry.yaml  # refresh the packaged default
```

Rules the registry layer enforces:

- Entries with missing, non-finite or out-of-range coordinates are skipped and
  listed. They are never set to zero (`bw-biust` publishes no altitude and is
  skipped). A duplicated TART id is ambiguous, so every copy is skipped.
  `--strict` makes any such entry fatal.
- `api_url` is the convention `https://api.elec.ac.nz/tart/<tart_id>` unless
  the registry publishes one. It does not mean the API is online.
- **Registry presence ≠ archived data.** `archive_name` is set only when
  `--verify-archives` found the prefix in the archive bucket, and
  `archive_verified_utc` records when. Verifications carry over from the
  previous catalogue, and a failed check never erases one.
- Each refresh is diffed against the previous catalogue: new/removed stations,
  name and metadata changes, and coordinate changes with their deltas and
  displacement in metres (WGS84 → Cartesian with Skyfield). A displacement
  above `--displacement-warning-m` (100 m by default) is logged as a warning.
  This is a flag level only. Nothing is rejected on that basis.
- If the registry cannot be fetched, `tart-stations-update` exits 1 and
  changes nothing.

**Reproducibility.** The live registry is for discovery only. `tart-common-gnss`
reads the frozen packaged catalogue by default and needs no network access.
`--refresh-stations` is opt-in, and it writes a run-specific catalogue under
`data/stations/`. Every run's metadata copies the selected stations' exact
coordinates and records the catalogue's registry snapshot hash. Stage 2 never
refreshes stations. It takes identity and coordinates from the Stage-1
results, and uses the catalogue only for archive and API metadata. If the
catalogue coordinates now differ, the manifest records
`station_metadata_changed: true` together with the `catalogue_coordinates`.

Packaged catalogue (registry retrieved 2026-10-09; 17 valid of 18 entries,
archive prefix verified for 14). Run `tart-stations-list` for the full list.
The seven stations below have the same coordinates as the 2026-10-07 legacy
file:

| TART id | Alias | Registry name | Lat (°) | Lon (°) | Alt (m) |
|---|---|---|---|---|---|
| `mu-udm` | `mauritius` | Bel Air - Mauritius | −20.2587508 | 57.7591989 | 20.0 |
| `mg-tana` | `madagascar` | Madagascar - Antananarivo | −18.8971446 | 47.5551677 | 1280.0 |
| `za-rhodes` | `south_africa` | Rhodes University, Grahamstown | −33.319413 | 26.507688 | 642.26 |
| `za-hammanskraal` | `south_africa_hammanskraal` | Kwalata - South Africa | −25.38948 | 28.32231 | 1163.0 |
| `na-unam` | `namibia` | Namibia - UNAM | −22.612053 | 17.856784 | 1655.0 |
| `zm-cbu` | `zambia` | Copperbelt University - Zambia | −12.80527778 | 28.23722222 | 1270.0 |
| `ghana` | `ghana` | GRAO - Ghana | 5.749725 | −0.305557 | 72.0 |

## Stage 1: common visibility (`tart-common-gnss`)

```bash
tart-common-gnss \
    --stations mauritius madagascar south_africa \
    --start 2026-10-01T00:00:00Z \
    --end   2026-10-02T00:00:00Z \
    --min-elevation 20 --step 10
```

`python common_gnss.py …` and `python -m tart_gnss …` are equivalent.

| Option | Meaning | Default |
|---|---|---|
| `--stations` | exactly three stations: canonical TART ids (`mu-udm`) or aliases (`mauritius`); see `tart-stations-list` | required |
| `--start`, `--end` | ISO-8601 times. `Z`/offsets are honoured; **naive times are interpreted as UTC** | required |
| `--min-elevation` | elevation threshold, degrees | 20 |
| `--step` | coarse sampling step, seconds | 10 |
| `--tolerance` | boundary refinement tolerance, seconds | 0.1 |
| `--tle-file` | use this TLE catalogue exactly | download from CelesTrak |
| `--tle-dir` | where downloaded catalogues are archived | `data/tle` |
| `--constellations` | restrict to e.g. `GPS Galileo` (`GPS`, `Galileo`, `BeiDou`, `GLONASS`) | all |
| `--output` | CSV path | `results/common_gnss_<YYYYMMDD>.csv` |
| `--station-file` | YAML/JSON station catalogue | packaged frozen registry snapshot |
| `--aliases` | station alias file | bundled `station_aliases.yaml` |
| `--refresh-stations` | fetch the live TART registry first and use the refreshed catalogue (archived under `--station-snapshot-dir`) | off |
| `--station-snapshot-dir` | where `--refresh-stations` writes snapshots | `data/stations` |
| `--min-duration` | drop windows shorter than this, seconds | 0 |
| `--tle-age-warning-days` | warn when \|t_obs − t_TLE epoch\| exceeds this (a flag level, not an accuracy bound) | 3 |
| `--verbose` | debug logging | off |

**Outputs** (shown for the default name):

| File | Contents |
|---|---|
| `results/common_gnss_<tag>.csv` | One row per satellite × window: rank, times, duration, clipping flags, and az/el at start/mid/end plus mid-window range for each station. Also the mid-window satellite position in **ITRS** (km) and the TLE epoch, age and source. |
| `results/common_gnss_<tag>.json` | Everything above plus full TLE lines, range rates, ITRS position and velocity at start/mid/end, line-of-sight unit vectors, station ITRS vectors and every searched satellite. **This is the Stage-2 input.** |
| `results/run_metadata_<tag>.json` | Software versions, command line, the selected stations (full coordinates), the station catalogue path and SHA-256, the registry source / retrieval time / snapshot SHA-256 (`station_registry_*`), thresholds, TLE file with SHA-256 and retrieval time, and a TLE-age summary. |

Results and metadata carry `"schema_version": 2`. Version 2 identifies
stations by canonical TART id (`mu-udm`) when the default catalogue is used;
the older unversioned files (e.g. those in `results/`) use catalogue names
such as `mauritius`, and Stage 2 still reads them.

**Method:** propagate each satellite on a coarse grid. Common visibility is
`min_i E_i(t) ≥ E_min`. Split the visible samples into maximal runs, refine
each rise and set by bisection on the exact predicate, then rank by duration.

## Archive scan: finding archived overlaps (`tart-archive-scan`)

Stage 1 finds *geometric* visibility, and a Stage-2 dry run checks one
window at a time. Before retrieving anything, check where the stations
actually archived data at the same time:

```bash
# 1. archive inventory only (no Stage-1 input)
tart-archive-scan --stations madagascar south_africa mauritius \
    --start 2026-10-09T00:00:00Z --end 2026-11-03T00:00:00Z \
    --min-common-duration 300 --output results/archive_candidates.csv

# 2. match every Stage-1 window, reusing the cached listings
tart-archive-scan --visibility-results results/common_gnss_20261009.json \
    --signal-relevance GPS_L1 GALILEO_E1 BEIDOU_B1C \
    --cache-policy reuse --output results/archive_candidates.csv
```

How it works:

- **Listing.** Each `(station, UTC day)` prefix is listed **once**, using
  the Stage-2 archive backend, prefix builder and filename parser. Only
  metadata is read; no file is ever downloaded.
- **Status per day.** `FOUND`, `EMPTY`, `ERROR` or `UNKNOWN`. A listing
  error never counts as an empty day.
- **Candidate coverage.** A file named `T` is *assumed* to cover
  `[T − 120 s, T]`. Gaps are never filled. The stations are intersected,
  and the shared segments are matched with all Stage-1 windows.
- **Labels.** Every row is labelled `ARCHIVE_COVERAGE_CANDIDATE` /
  `FILENAME_CANDIDATE`, and `verified_overlap_duration_s` stays null.
- **Ranking.** Rows are ranked by:
  1. complete inventory;
  2. in-band signal;
  3. duration, capped at 600 s;
  4. |TLE age|;
  5. Stage-1 geometry;
  6. Stage-1 rank.
- **Handoff.** Each row gives `tart-retrieve --visibility-results … --rank N --dry-run`.

**Exit status:**

- `0`: complete inventory, even with zero candidates;
- `3`: inventory incomplete (listing errors or unknown days);
- `2`: invalid input or a station set that differs from Stage 1;
- `1`: a failure, or every listing failed.

**Outputs:** `archive_candidates.csv/.json`, `archive_candidate_stations.csv`,
`archive_common_segments.csv`, `archive_station_days.csv`,
`archive_station_summary.csv` and `archive_scan_manifest.json`. The cache is
written to `data/archive_inventory/` (git-ignored).

**2026-10-09 → 11-03 result:** `mu-udm` has no archived files at all (its
archive ends on 2026-06-26), so with a complete inventory there is **no**
three-station overlap. A `--min-stations 2` diagnostic finds 131 h of
mg-tana + za-rhodes candidate time. Method, cache, ranking, real results
and the decision log are in [`docs/archive_scan.md`](docs/archive_scan.md).

### Network-wide discovery (`--all-stations`)

Instead of guessing a triad, scan the whole registry once:

```bash
tart-archive-scan --all-stations \
    --start 2026-07-01T00:00:00Z --end 2026-11-03T00:00:00Z \
    --min-stations 3 --min-common-duration 300 --min-baseline-km 500 \
    --cache-policy reuse --output results/network_overlap.csv
```

- **Stations** come from the registry at run time: an `archive_name` and
  valid coordinates are required. Excluded entries are reported with a
  reason (`NO_ARCHIVE_NAME`, `INVALID_COORDINATES`, ...).
- **Physical sites.** Stations within `--site-tolerance-km` (default 1 km)
  form one site, so co-located receivers (e.g. `nz-nelson`/`nz-nelson-32`)
  never count as independent. `--min-stations` counts **sites**.
- **One sweep** over all stations' candidate segments gives the network
  intervals with at least N sites. Triads (stations at distinct sites) are
  then built only inside those intervals, each from its *own* members'
  intersected coverage, so a triad stays one segment while a fourth station
  joins or leaves.
- **Baselines** are WGS84 geodesics (Vincenty); `--min-baseline-km` requires
  *every* pair to pass. It is a site-diversity filter, not a triangulation
  metric.
- **Ranking:** complete inventory, total possible shared time, longest
  contiguous segment, minimum baseline, station ids.
- **Outputs:** `network_overlap.csv/.json` (ranked triads),
  `network_active_segments.csv`, `network_triad_segments.csv`,
  `network_station_days.csv`, `network_station_summary.csv`,
  `network_excluded_stations.csv` and `network_archive_scan_manifest.json`
  (with output SHA-256s).
- **Handoff.** Each triad gets a *proposed* `tart-common-gnss` command for
  its longest segment and the follow-up explicit-station scan. Nothing is
  run, and `--visibility-results` is rejected in this mode, because a
  Stage-1 JSON belongs to one station set.

**2026-07-01 → 11-03 result.**

- **Stations.** 14 eligible stations at 13 sites, with all 1,778
  station/day listings complete. 5 of the 14 stations have no files in the
  interval (`mu-udm`, `na-unam`, `nz-dunedin`, `se-kth`, `zm-cbu`).
- **Candidates.** 44 triads have candidate overlap, and all of them pass the
  500 km cut.
  - Top triad: `bd-iub` + `ghana` + `za-rhodes`, with 1,750 h of possible
    shared time and a longest segment of 91 h.
  - `bd-iub` + `mg-tana` + `za-rhodes` ranks 5th, with 760 h.
- These are archive candidates only. Stage 1 must still establish common
  GNSS visibility. Details are in
  [`docs/archive_scan.md` §9](docs/archive_scan.md#9-network-wide-discovery---all-stations).

## Stage 2: observation retrieval (`tart-retrieve`)

Check what the archive holds before downloading anything:

```bash
tart-retrieve --visibility-results results/common_gnss_20260625.json --rank 8 --dry-run
```

Retrieve a short subsection first, then a longer interval:

```bash
tart-retrieve \
    --visibility-results results/common_gnss_20260625.json \
    --rank 8 --padding 120 --max-duration 600 \
    --output-dir data/observations \
    --download --convert-ms --keep-native
```

| Option | Meaning | Default |
|---|---|---|
| `--visibility-results` | Stage-1 JSON. A Stage-1 CSV is accepted if its JSON sits next to it | required |
| `--rank` / `--norad-id` / `--window-index` | select exactly one window. A NORAD ID with several windows needs `--window-index` (1-based, time order) | — |
| `--padding` | seconds added before and after the predicted window | 120 |
| `--max-duration` | retrieve at most this many seconds, centred on the window midpoint | whole window |
| `--output-dir` | root output directory | `data/observations` |
| `--download` | download the native HDF5 files | off |
| `--convert-ms` | convert each station to its own Measurement Set (needs `--download`) | off |
| `--keep-native` | keep native files after a validated conversion (their hashes are always kept) | off |
| `--allow-partial` | continue when a station has no data. The result **stays invalid** | off |
| `--uncalibrated` | ask tart2ms not to apply the gains/phases stored in the files | off |
| `--tart2ms` | path to the tart2ms executable | auto-detect |
| `--gap-factor` | a step longer than this × the median cadence counts as a gap | 3 |
| `--dry-run` | select, resolve and query only; write nothing | off |
| `--station-file`, `--aliases`, `--verbose` | as in Stage 1 (the catalogue supplies only archive/API metadata) | — |

**Exit status:** `0` means a validated three-station dataset; `3` means the
run completed but is **not** three-station valid (the reasons are printed
and stored in the manifest); `2` means invalid input or selection; `1` means
a failure.

**Output layout:**

```text
data/observations/20260625/NORAD_40534/     # _w<k> is appended if the satellite has several windows
├── stage1_window.json        # the selected Stage-1 window, unchanged
├── manifest.json             # Stage-2 → Stage-3 interface
├── retrieval.log
├── mauritius/raw/vis_*.hdf (+ .prov.json)   mauritius/ms/mauritius.ms (+ tart2ms logs)
├── madagascar/raw/…                          madagascar/ms/madagascar.ms
└── south_africa/raw/…                        south_africa/ms/south_africa.ms
```

**How it works**

- **Archive.** The public TART S3 bucket (`tart`, `s3.us-west-2.amazonaws.com`;
  location reused from `tart_tools.archive_handler`), laid out as
  `<archive_name>/vis/<Y>/<M>/<D>/vis_<UTC>.hdf`. Files are selected by the
  observation time in the filename, which is the last-sample time. They are
  *not* selected by S3 upload time, which lags by minutes.
- **Integrity.** Original filenames are kept. Size and the MD5/ETag are
  checked, and the SHA-256 goes into a `.prov.json` sidecar. A verified local
  copy is reused; any other existing file stops the run.
- **Validation.** Each file is read: timestamps, `vis` shape vs baselines,
  antennas (24) and baselines (276), and frequency and bandwidth from the
  file's own `config`. The frequency counts as in band if 1575.42 MHz lies
  within the recorded band. The station position must match the
  Stage-1 station coordinates, and zero-gain antennas and calibration changes are reported.
- **Synchronisation.** The real common interval is `[max first, min last]`,
  and it must have positive length. Coverage fractions of the predicted
  (and retrieved) window count gaps as uncovered. Nearest-sample offsets are
  given for each station pair; the stations are not assumed to sample
  simultaneously.
- **Measurement Sets.** `tart2ms --hdf … --ms … --skip-online-source-catalog
  --skip-celestial-source-catalog` runs once per station, with an argument
  list rather than a shell string. Each MS is validated with python-casacore:
  rows, TIME, UVW, DATA, ANTENNA, SPECTRAL_WINDOW. Stations are **never**
  combined into one interferometer.
- **Calibration.** By default tart2ms applies the gains and phases stored
  in each file and flags zero-gain antennas. The manifest records this as
  *applied, unverified*.

**`manifest.json`** has the sections `three_station_valid` and
`validity_reasons`, `satellite` (including `signal_relevance`:
`GPS_L1`/`GALILEO_E1`/`BEIDOU_B1C`/`OUT_OF_BAND`/`UNKNOWN`), `prediction`,
`stations`, `synchronisation`, `stage1_window`, `software` and `provenance`.
Each `stations` entry holds the archive status, the raw files with hashes,
the native metadata, coverage, calibration, MS metadata and the predicted
az/el. A previous manifest is renamed, never overwritten.

## Stage 3: direction estimation (`tart-localise`)

```bash
pip install -e ".[stage3]"      # scipy, pandas, pyarrow, matplotlib, python-casacore
tart-localise \
    --manifest data/observations/20260625/NORAD_40534/manifest.json \
    --all-stations --method visibility --image-validation --time-bin 1 \
    --output-dir results/stage3
```

| Option | Meaning | Default |
|---|---|---|
| `--manifest` | Stage-2 `manifest.json` (authoritative input) | required |
| `--station` / `--all-stations` | one or more stations (`station_id` or `tart_id`), or all, each processed independently | — |
| `--method` | primary estimator: `visibility`, `image`, or `both` | `visibility` |
| `--image-validation` | also run the dirty-image estimator as a cross-check | off |
| `--time-bin` | bin width (s); empty bins are skipped, never fabricated | 1 |
| `--bin-origin` | `utc`: edges on a common UTC grid (same bins at every station); `first-sample`: from each station's first sample | `utc` |
| `--min-baselines` | fewer valid baselines → no estimate, `LOW_BASELINE_COUNT` | 10 |
| `--search-radius-deg` | search cap about the TLE direction (search bound only, never a prior) | 10 |
| `--blind` | search the whole visible hemisphere | off |
| `--robust-loss` | `linear`, `huber`, `soft_l1` | `huber` |
| `--bootstrap N`, `--seed` | baseline bootstrap resamples per bin and its seed | 0, 20260625 |
| `--min-snr`, `--max-flagged-fraction` | quality thresholds | 5, 0.5 |
| `--association-radius-deg` | post-fit check: fits further than this from the TLE are flagged `TLE_UNASSOCIATED` (≤ 0 disables) | 5 |
| `--allow-invalid-stage2` | process valid stations of a non-three-station-valid run (marked `stage2_three_station_valid: false`) | off |
| `--data-root` | directory the manifest's relative MS paths refer to | auto |
| `--no-plots`, `--verbose` | skip figures; debug logging | off |

**Exit status:** `0` means at least one usable estimate; `3` means it completed without any; `2` means invalid input; `1` means a failure.

**Output:** `results/stage3/<date>/NORAD_<id>/` with `directions.csv`,
`directions.parquet`, `stage3_manifest.json`, `diagnostics/*.png` and
`logs/stage3.log`. Each row has the fitted (l, m) and (Az, El) with
covariances, complex amplitude, χ², residual statistics, SNR, baseline
counts, the image-domain estimate and separation, optional bootstrap
percentiles, `usable` and `quality_flags`. TLE az/el columns are reference
predictions only.

**Conventions.** ENU frame, Az from North through East, l = East and m = North
direction cosines. The model is V = C exp(−2πi ν/c (r<sub>q</sub> − r<sub>p</sub>)·s). tart2ms stores
UVW = r<sub>p</sub> − r<sub>q</sub> for a zenith phase centre. This sign was validated on the real
data (a mirrored convention would put PRN 26 and PRN 28 at Az + 180°), and
every MS is checked against it on load. The near-field metric B²/(λR) ≈ 10⁻⁵
is recorded per station, so the plane-wave model holds. Full derivation,
uncertainty method, flags and limitations:
[`docs/stage3_direction_estimation.md`](docs/stage3_direction_estimation.md).

## Stage 4: 3-D localisation (`tart-triangulate`)

```bash
tart-triangulate \
    --stage3-manifest results/stage3/20260625/NORAD_40534/stage3_manifest.json \
    --output-dir results/stage4
```

| Option | Meaning | Default |
|---|---|---|
| `--stage3-manifest` | `stage3_manifest.json` from `tart-localise` | required |
| `--time-match-tolerance` | max \|t<sub>station</sub> − t<sub>reference</sub>\| (s); nearest row, no interpolation | 0.75 |
| `--min-stations` | stations needed for a valid estimate (2-station solutions are degraded diagnostics) | 3 |
| `--method` | `weighted` (Stage-3 covariance × range²) or `unweighted` (closed form Σ P<sub>i</sub> x = Σ P<sub>i</sub> r<sub>i</sub>) | `weighted` |
| `--quality-policy` | `strict`: only Stage-3 `usable` rows; `permissive`: also `BOUNDARY_SOLUTION`/`LOW_SNR`/`TLE_UNASSOCIATED` rows (flagged) | `strict` |
| `--bootstrap N`, `--seed` | parametric bootstrap of the directions per epoch | 0, 20260625 |
| `--stations`, `--reference-station`, `--reference-time` | restrict stations; choose the reference clock; one epoch only | all, most rows, — |
| `--max-condition-number`, `--max-normalised-residual`, `--max-reduced-chi2`, `--time-mismatch-sigma` | flag thresholds | 1e8, 5, 5, 1 |

**Exit status:** `0` at least one valid 3-D estimate; `3` none; `2` invalid input; `1` failure.

**Output:** `results/stage4/<date>/NORAD_<id>/` with `positions.{parquet,csv}`
(x, y, z, the six covariance terms, σ, station count, line residuals, χ²,
condition number, crossing angles, timing, bootstrap percentiles, geodetic
sub-point, flags), `station_contributions.{parquet,csv}` (per station and epoch:
time offset, ECEF station position, Az/El and covariance, range, line
residual, angular rate and motion during the offset, light time, Stage-3 and
Stage-4 flags), `stage4_manifest.json`, five diagnostic plots and
`logs/stage4.log`. Station coordinates come from the run's frozen Stage-1/2
provenance (cross-checked against the hash-verified Stage-2 manifest when it
is available), never from the live registry.

The linear position covariance agrees with Monte Carlo scatter to within a
few per cent for angular errors up to about 0.5°. It is long and thin along
the line of sight: with 1 000–3 000 km baselines at about 20 000 km range, the
along-LOS σ is 5–10× the transverse σ. Method, flags and validation:
[`docs/stage4_triangulation.md`](docs/stage4_triangulation.md).

## Stage 5: trajectory tracking (`tart-track`)

```bash
pip install -e ".[stage5]"      # scipy, pandas, pyarrow, matplotlib
tart-track \
    --stage4-manifest results/stage4/20260625/NORAD_40534/stage4_manifest.json \
    --dynamics two-body \
    --output-dir results/stage5
```

| Option | Meaning | Default |
|---|---|---|
| `--stage4-manifest` | `stage4_manifest.json` from `tart-triangulate` | required |
| `--dynamics` | `two-body` (`solve_ivp` + variational equations), `constant-velocity` (diagnostic baseline) or `j2` | `two-body` |
| `--process-accel-sigma` | white process acceleration σ (m/s²) | per model: 1 (CV), 1e-4 (two-body), 1e-5 (J2) |
| `--initial-velocity-window` | seconds of Stage-4 data used for the initial (r, v) fit | 30 |
| `--innovation-gate-probability` | χ²₃ probability of the innovation gate | 0.9973 |
| `--quality-policy` | `strict`: Stage-4 `valid` epochs not built on unusable Stage-3 rows; `permissive`: also `POOR_GEOMETRY` and `STAGE3_QUALITY_WARNING` epochs | `strict` |
| `--max-gap-s` | `LARGE_TIME_GAP` after this long without an assimilated measurement | 120 |
| `--robust-update` | `none`: gated measurements are not assimilated; `inflate`: they are, with covariance inflated to the gate | `none` |
| `--no-smoother`, `--no-plots`, `--verbose` | filtered states only; skip figures; debug logging | off |

**Exit status:** `0` means ≥ 3 measurements assimilated and ≥ 1 valid row;
`3` means it completed without enough valid estimates; `2` means invalid input;
`1` means a failure.

**Output:** `results/stage5/<date>/NORAD_<id>/` contains:

- `trajectory.{parquet,csv}`, one row per Stage-4 epoch: filtered and
  smoothed GCRS position and velocity, the smoothed ECEF position and its
  3×3 covariance, speed, σ, model acceleration, the Stage-4 measurement,
  used/accepted, NIS, `valid` and flags, plus the full filtered and smoothed
  6×6 covariances as 21 upper-triangle columns each;
- `innovations.{parquet,csv}`: ν, S, whitened ν, NIS, gate and accepted;
- `stage5_manifest.json`;
- seven plots in `diagnostics/`;
- `logs/stage5.log`.

**How it avoids false precision:**

- The Stage-4 covariance is kept in full.
- Process noise is explicit and recorded.
- Gaps widen the prediction.
- Gated measurements stay in the output.
- The initial fit is used as a 10⁴-inflated prior, so the initial window is
  not counted twice.
- Covariances are checked after every step. Invalid Stage-4 covariances are
  rejected, never repaired.

After 3 consecutive gated measurements the filter restarts from the data
(`FILTER_DIVERGENCE`) instead of locking out.

On synthetic orbits the predicted state σ matches the Monte Carlo error to
within about 2%, and the coverage is close to χ²₃. Method, flags and validation
are in [`docs/stage5_tracking.md`](docs/stage5_tracking.md).

## Stage 6: TLE validation (`tart-validate`)

```bash
pip install -e ".[stage6]"      # pandas, pyarrow, matplotlib (SGP4 via skyfield/sgp4)
tart-validate \
    --stage5-manifest results/stage5/20260625/NORAD_40534/stage5_manifest.json \
    --tle-file data/tle/gnss_20260625.tle \
    --output-dir results/stage6
```

| Option | Meaning | Default |
|---|---|---|
| `--stage5-manifest` | `stage5_manifest.json` from `tart-track` | required |
| `--tle-file` / `--tle-source` | exactly one primary reference: a TLE file (nearest epoch if it holds several sets for the NORAD ID), `stage1` (the exact Stage-1 set, from the run provenance) or `download` (explicit CelesTrak download, archived; current elements only) | required |
| `--tle-ensemble-files`, `--tle-ensemble-dir`, `--tle-ensemble-window-days` | nearby sets for the *empirical reference dispersion* (≥ 3 needed) | none, 7 d |
| `--tle-age-warning-days`, `--tle-age-stale-days` | TLE-age thresholds (pipeline quality limits, not SGP4 accuracy) | 3, 14 |
| `--quality-policy` | `strict`: valid Stage-5 rows; `permissive`: any row with a finite state | `strict` |
| `--position-residual-warning-km`, `--velocity-residual-warning-m-s` | optional residual flags (no built-in thresholds) | off |
| `--allow-stale-reference` | exit 0 for a descriptive comparison when TLE age is the only failed condition | off |
| `--no-plots`, `--verbose` | skip figures; debug logging | off |

**Exit status:** `0` means precision validation passed; `3` means the
comparison completed but is not valid for precision analysis (for example a
stale TLE); `2` means invalid input, including `TLE_NOT_FOUND` and
`TLE_NORAD_MISMATCH`; `1` means a failure, including `TLE_FRAME_FAILURE`.

**Output:** `results/stage6/<date>/NORAD_<id>/` contains:

- `residuals.{parquet,csv}`: observed and TLE states, Δr and Δv (GCRS and
  RTN), RTN σ of the observation (and of the reference with an ensemble),
  the angular and range residual, the along-track time equivalent, squared
  Mahalanobis distances, TLE age and flags;
- `tle_trajectory.parquet`: the TLE states and RTN basis vectors, so Stage 7
  need not re-propagate;
- `summary.json`: mean/median/RMS/MAD/percentiles, trends, autocorrelation
  and the precision-validation conditions;
- `stage6_manifest.json`;
- eight plots;
- `logs/stage6.log`.

The TLE has no formal covariance. `*_obs_cov` distances are therefore
conditional on the TLE being fixed, and reference-aware distances need an
ensemble of nearby TLEs (labelled empirical dispersion).

The frame chain (SGP4 TEME → Skyfield → GCRS) agrees with an independent
route (raw SGP4 → GMST82 → the Stage-5 rotation) to 2.5 cm. Method, flags and
validation are in [`docs/stage6_validation.md`](docs/stage6_validation.md).

## Stage 7: network diagnosis (`tart-diagnose`)

```bash
pip install -e ".[stage7]"      # scipy, pandas, pyarrow, matplotlib
tart-diagnose \
    --stage6-manifest results/stage6/20260625/NORAD_40534/stage6_manifest.json \
    --output-dir results/stage7
```

| Option | Meaning | Default |
|---|---|---|
| `--stage6-manifest` | `stage6_manifest.json` from `tart-validate`; Stages 5–2 are found via the recorded provenance and verified by SHA-256 | required |
| `--quality-policy` | `strict`: only valid, fully linked epochs, and hash mismatches are fatal; `permissive`: also rows with a finite residual, mismatches recorded | `strict` |
| `--leave-one-out` / `--no-leave-one-out`, `--pairwise` / `--no-pairwise` | re-triangulation analyses | on |
| `--correlation-method` | `spearman`, `pearson` or `kendall` | `spearman` |
| `--min-samples`, `--min-correlation` | sample-count and effect-size requirements for a correlation to support a tag | 8, 0.5 |
| `--anomaly-threshold` | Stage-6 squared-Mahalanobis anomaly threshold | χ²₃ at 1 − 0.0027/n_epochs |
| `--min-persistent-epochs` | consecutive strong epochs needed for persistence (Level 5) | 5 |
| `--alignment-tolerance-s`, `--no-plots`, `--verbose` | time-join tolerance; skip figures; debug logging | 0.01 s |

**Evidence levels** describe evidence patterns, never causes:

| Level | Meaning |
|---|---|
| 0 | normal |
| 1 | statistical anomaly |
| 2 | station-specific |
| 3 | instrument, calibration, timing or localisation correlated |
| 4 | network-consistent: no station subset removes it, and it is not geometry-, trajectory-model- or reference-limited |
| 5 | unexplained: Level 4, strong and persistent on the epoch-by-epoch Stage-4 data, with complete diagnostics and a known TLE reference uncertainty (ensemble) |

A stale reference caps the level at 3. Missing calibration metadata, or a
single TLE with unknown uncertainty, caps it at 4. Confidence is a transparent LOW / MODERATE / HIGH from six
components, never a probability.

**Exit status:** `0` means the diagnosis completed; `3` means fewer than 3
diagnosable epochs; `2` means invalid input or a broken provenance chain;
`1` means a failure.

Method, rules, validation and the real-data run are in
[`docs/stage7_diagnosis.md`](docs/stage7_diagnosis.md).

## Campaign analysis (`tart-campaign`)

```bash
pip install -e ".[campaign]"    # scipy, pandas, pyarrow, matplotlib
tart-campaign \
    --input-root results/stage7 \
    --output-dir results/campaign \
    --campaign-name gnss_validation_2026
```

| Option | Meaning | Default |
|---|---|---|
| `--input-root` / `--stage7-manifest` | recursive search for `stage7_manifest.json` and/or explicit manifests (repeatable; at least one required) | — |
| `--quality-policy` | `strict`: runs with broken provenance are excluded; `permissive`: kept and flagged | `strict` |
| `--precision-only` / `--no-precision-only` | TLE-relative summaries from the precision-validation cohort only | on |
| `--precision-validation-only` | exclude non-precision runs from the campaign | off |
| `--min-epochs-per-run`, `--min-runs-per-group` | run and group size requirements | 3, 3 |
| `--group-by` | extra grouped tables: `constellation`, `norad_id`, `station_set`, `cohort`, `observation_date`, `tle_source` | `constellation` |
| `--reference-frame` | residual components summarised: `rtn` or `gcrs` | `rtn` |
| `--bootstrap N`, `--seed` | run-level bootstrap resamples and seed | 1000, 20260625 |
| `--norad-id`, `--constellation`, `--station`, `--start`, `--end`, `--min-elevation`, `--max-tle-age-days`, `--max-evidence-level` | filters; every excluded run is listed with its reason | — |
| `--publication-plots`, `--no-plots`, `--verbose` | larger fonts and vector PDFs; skip figures; debug logging | off |

**Exit status:** `0` means the analysis completed; `3` means fewer than 3
runs in the analysis cohort (descriptive outputs are still written,
`INSUFFICIENT_CAMPAIGN_SIZE`); `2` means invalid input; `1` means a failure.

**Output:** `results/campaign/<name>/` contains `campaign_runs`,
`campaign_epochs` and `campaign_stations` (Parquet + CSV), `excluded_runs.csv`,
`summary.json`, `campaign_manifest.json`, `tables/*.csv` and ten figures.

Residuals are reported as **TLE-relative**, never as absolute accuracy, and
no physical-cause rates are produced. Method, statistics and validation are
in [`docs/campaign_analysis.md`](docs/campaign_analysis.md).

## Simulation validation (tab-sim, research)

A research track, kept separate from the pipeline, checks direction recovery
against visibilities simulated with
[tab-sim](https://github.com/chrisfinlay/tab-sim) 0.0.13: a single GNSS
satellite (NORAD 39166, L1) at finite range, observed by the Ghana TART
(24 antennas, 120 × 1.024 s). The fit is orbit-blind: a spherical-wave,
antenna-*w*-corrected forward model with a profiled complex amplitude.
The saved orbit is used only to score the result. **Stages 1–7 are not
changed**, and nothing here is an accuracy claim for real TART data.

| Check | Result |
|---|---|
| Forward-model phases (orbit truth) | Finite range + antenna *w*: coherence ≈ 1.000000, circular RMS 0.048°. Without *w*: coherence 0.038. |
| Noise-free orbit-blind fit | Az 100.9484°, El 51.7286°; residual 5.8 × 10⁻¹²; about 3.7″ from the OMM direction. Insensitive to the assumed range (18 000–26 000 km). |
| The 3.7″ offset | Traced to time scales: tab-sim passes the UTC-labelled Julian dates to Skyfield's `ts.ut1_jd`. Against a simulator-consistent (UT1) reference the offset is (−3.548″, +1.097″) and the profile residual drops from 3.8 × 10⁻⁷ to 1.2 × 10⁻¹². It is a truth-reference convention, not an interferometric bias. The upstream simulator has not been fixed. |
| Monte Carlo, 200 trials per SNR (1 200 fits), per-visibility SNR 3–100 | Median error against the simulator reference: 111″ (SNR 3), 30″ (SNR 10), 3.4″ (SNR 100). At SNR 100 the bias falls from (−3.57″, +0.88″) to (−0.03″, −0.22″). Raw held-out 68 % ellipse coverage is 0.54–0.74, so the empirical covariance is not yet a calibrated posterior. |

Commands, environment and outputs are in
[`docs/tabsim_gnss/README.md`](docs/tabsim_gnss/README.md). The report is
[`docs/report/tabsim_gnss_validation.pdf`](docs/report/tabsim_gnss_validation.pdf)
(source `tabsim_gnss_validation.tex`). The scripts are in `scripts/tabsim_*.py`, and the simulation inputs and
run records (per-station tab-sim YAML, ITRF antenna positions, resolved
configs, the OMM record used and logs) are in [`configs/tart_tabsim/`](configs/tart_tabsim).
They need astropy, scipy, skyfield and python-casacore, plus tab-sim to
regenerate the simulations. The Measurement Sets and results are local and
not committed.

## TLEs and reproducibility

- Without `--tle-file`, current GNSS TLEs are downloaded from CelesTrak and
  archived as `data/tle/gnss_<retrieval>.tle`, with a `.meta.json` sidecar
  holding the source URLs, retrieval time, SHA-256 and constellations.
  Archived catalogues are never overwritten, and a modified file is rejected.
- **CelesTrak serves only current elements.** For older observations, supply
  a catalogue from the observation date with `--tle-file`. Stage 1 warns when
  it propagates backwards or far from a TLE epoch.
- Stage-2 manifests record the SHA-256 of the Stage-1 results and of the
  station file, the archive identifiers and SHA-256 of every native file, the
  exact tart2ms command and all software versions.

## Example results

| Run | Result |
|---|---|
| Stage 1, 2026-10-01 ([`results/`](results)) | 154 common windows, 135 of 145 satellites. |
| Stage 2, spec acceptance command (2026-10-01, rank 1) ([example](docs/examples/stage2/invalid_20261001_NORAD_38091)) | **Not valid** (exit 3). BeiDou-2 G5 is `OUT_OF_BAND`, and **Mauritius has no archive data on that day**: the `mu-udm` archive ends on 2026-06-26. Nothing was downloaded. |
| Stage 1, 2026-06-25, GPS+Galileo ([`results/`](results)) | 70 windows. All three stations have 24-h archive coverage that day. |
| Stage 3, 2026-06-25 PRN 26, real data, 10 s bins ([details](docs/stage3_direction_estimation.md#17-real-data-run-2026-06-25-gps-prn-26-norad-40534)) | ZA-Rhodes tracks PRN 26 in all 68 bins (SNR ≈ 12), within 1° of the 103-day-old TLE. MU-UDM is usable from 09:18:36 onwards (37 bins); before that a brighter source pulls the fit to the cap edge (flagged). MG-Tana is dominated by a quasi-static non-target source (flagged `TLE_UNASSOCIATED`). The bin-to-bin scatter (0.02–0.05°) is within the Jacobian σ (0.06–0.15°). The real error is a slowly varying systematic (0.1–0.2° wander, 1–2° offset from the TLE), which the bootstrap (0.4–0.7°) partly reflects. |
| Stage 4, synthetic end-to-end ([example](docs/examples/stage4)) | A 3-D truth placed 1000 km off the PRN 26 TLE track (across the line of sight) is observed via simulated visibilities → `tart-localise` → `tart-triangulate`. The recovered position lies a median 996 km from the TLE along the injected offset (scatter 35 km). The 3-D error has median 293 km, almost all along the line of sight (transverse 32 km), and the pulls have std 0.87–0.99. |
| Stage 4, 2026-06-25 PRN 26, real data, 10 s UTC-aligned Stage-3 bins ([details](docs/stage4_triangulation.md#12-real-data-run-2026-06-25-gps-prn-26-norad-40534), [example](docs/examples/stage4_real)) | **The spec command exits 0 with 1 valid three-station epoch** (09:23:54.9). Stations now match within 0.57 s; MG-Tana's 2 usable rows are the limit. With two stations allowed there are 35 epochs: crossing angles 4–8°, condition number about 240, σ ≈ (266, 163, 453) km. The solutions lie about 3 000 km too far out along the line of sight (height about 23 300 km), with reduced χ² ≈ 11, which matches the 1–2° Stage-3 direction systematics. This is not a TLE comparison. |
| Stage 5, synthetic ([example](docs/examples/stage5), [report](docs/examples/stage5/demo_report.txt)) | Circular GNSS-like orbit (analytic Kepler truth, no TLE). From 156 Stage-4-like positions (30 km along the line of sight, 3 km across, with a 240 s gap and a 150 km outlier), the smoothed trajectory has a position RMSE of 5.2 km and a velocity RMSE of 1.1 m/s, and the outlier is gated. Over 200 Monte Carlo runs, the predicted σ matches the empirical error (3.45 vs 3.53 km; 3.61 vs 3.62 m/s), with coverage close to χ²₃. Two-body beats constant velocity on the innovation diagnostics. End to end through `tart-triangulate`, a 39 km single-epoch error becomes a 6.9 km trajectory error. |
| Stage 5, 2026-06-25 PRN 26, real data ([details](docs/stage5_tracking.md#14-real-data-run-2026-06-25-gps-prn-26-norad-40534), [example](docs/examples/stage5_real)) | **The spec command exits 3** (`INSUFFICIENT_MEASUREMENTS`), because the default Stage-4 run has one valid epoch. On the two-station Stage-4 run (35 epochs over 340 s), the epoch-to-epoch scatter is consistent with the Stage-4 covariance (mean NIS 2.8, close to χ²₃). NIS cannot see a common offset, such as the ~3 000 km along-LOS offset found in Stage 4. The innovations are strongly autocorrelated (lag-1 0.69, `DYNAMICS_MISMATCH`), which shows slowly varying Stage-3 systematics, so the smoothed σ is optimistic. The only three-station epoch is gated (NIS 36). Two-body and constant velocity cannot be told apart at this noise level. No TLE comparison is made (that is Stage 6). |
| Stage 6, synthetic ([example](docs/examples/stage6), [report](docs/examples/stage6/demo_report.txt)) | With a fresh reference TLE (0.11 d), injected RTN offsets of (+1000, +2000, −500) m and Δv of (+0.10, −0.05, +0.02) m/s are recovered exactly, and the conditional D² equals its analytic value (18.92). End to end through `tart-track`, a 2 s time shift (7.818 km along track) is recovered as +7.789 km (σ 51 m), 1.992 s. With the TLE equal to the truth, the conditional position D² has mean 2.82 and median 2.33 over 100 Monte Carlo runs (χ²₃: 3, 2.37). |
| Stage 6, 2026-06-25 PRN 26, real data ([details](docs/stage6_validation.md#12-real-data-run-2026-06-25-gps-prn-26-norad-40534), [example](docs/examples/stage6_real)) | **Not a precision validation (exit 3).** The only local TLE (the Stage-1 set, epoch 2026-10-06) is propagated 103 days backwards: `TLE_STALE`, `TLE_BACKPROPAGATED`, `precision_validation: false`. The descriptive residuals (observed − TLE) are R +2765 km, T +370 km, N −627 km (median), about 1.4° geocentric. They are not interpreted as accuracy. A contemporaneous TLE (e.g. from Space-Track) is needed. |
| Stage 7, synthetic ([example](docs/examples/stage7), [report](docs/examples/stage7/demo_report.txt)) | **10 of 10 scenarios match their expected diagnosis**, each built through the real Stage 4–6 chain. A clean network gives Level 0 NORMAL. A 0.036° bias at one station is identified as Level 2 `STATION_SPECIFIC` (Madagascar, 36/36 epochs). A calibration-gated bias gives `CALIBRATION_CORRELATED`, and a late station gives `TIMING_CORRELATED`. Alternating geometry gives `GEOMETRY_CORRELATED` (Level 0), unmodelled acceleration gives `TRAJECTORY_MODEL_CORRELATED`, and a stale TLE gives `REFERENCE_LIMITED` (capped). A network-wide 10 s offset is Level 5 `NETWORK_CONSISTENT` + `UNEXPLAINED` with a TLE ensemble, and Level 4 with a single TLE or without Stage-2 metadata. Of 20 clean networks with random-walk metadata, none gives a spurious correlation tag. |
| Stage 7, 2026-06-25 PRN 26, real data ([details](docs/stage7_diagnosis.md#10-real-data-run-2026-06-25-gps-prn-26-norad-40534), [example](docs/examples/stage7_real)) | The whole chain to the Stage-2 manifest verifies, and re-triangulation reproduces Stage 4 to 0 m. The result is **Level 3, `REFERENCE_LIMITED`** (stale TLE), with `LOCALISATION_CORRELATED` (11 of 37 station pairs inconsistent) and `TRAJECTORY_MODEL_CORRELATED` (Stage-5 `DYNAMICS_MISMATCH`). Only one three-station epoch exists, so the influential station is undetermined. None of 130 correlation tests is supported. Madagascar's single direction is 1.4° off the network trajectory. |
| Campaign, synthetic ([example](docs/examples/campaign), [report](docs/examples/campaign/demo_report.txt)) | 17 runs built through the real Stage 3–7 chain: 3 satellites in 3 constellations, repeated passes, 2–4 stations and TLE ages 0.35–20 d. Stage-4 σ is 16.6 km (major/minor axis 16.5/1.1 km); Stage-5 σ is 4.5 km; the TLE-relative residual is 13.8 km (95 % CI 8.1–18.4, with injected age drift). Every injected property is recovered: the unstable station (Madagascar, only one tagged `STATION_SPECIFIC`, 18 % of runs), the geometry effect (σ 22.9/16.6/12.2 km for 2/3/4 stations), TLE-age drift (+1.13 km/day against 1.55 injected), the Galileo–GPS SNR difference (exact permutation p = 0.0018) and pass-to-pass scatter (3.58 against 3.87 km). |
| Campaign, real ([example](docs/examples/campaign_real)) | **Exit 3, `INSUFFICIENT_CAMPAIGN_SIZE`.** 2 runs were discovered: 1 is included (descriptive only, stale TLE) and 1 is excluded (`INSUFFICIENT_EPOCHS`). No population result is claimed. |
| Stage 3, synthetic end-to-end ([example](docs/examples/stage3)) | Sources injected 3.0–3.7° from the real PRN 26 TLE track are recovered to a median of 0.11–0.17° (σ ≈ 0.15°) and stay 3.0–3.7° from the TLE (measurement-driven). The visibility–image disagreement is about 0.3°. |
| Stage 2, 2026-06-25, rank 8 GPS BIIF-9, 10-min subsection ([example](docs/examples/stage2/valid_20260625_NORAD_40534)) | **Valid** (exit 0). Three validated 24-antenna MSs (182 160 rows each). Real overlap 09:13:39–09:24:24 (645 s, 100 % of the retrieved interval). Station dumps are offset by 0.21–0.46 s. **9 of 24 Rhodes antennas had zero gains** (62 % flagged), and Mauritius and Madagascar each changed gain solution once. |

These examples were produced before the commands were renamed, so their
recorded command lines use the old names `common-gnss` (in
`results/run_metadata_*.json`) and `retrieve-tart` (in the example
manifests). They are kept unchanged as provenance records.

## Caveats

- **Geometric visibility ≠ detectability.** GPS L1, Galileo E1 and BeiDou-3
  B1C are in the TART band. BeiDou-2 (1561 MHz) and legacy GLONASS FDMA
  (1598–1606 MHz) are not. Use `--constellations GPS Galileo` in Stage 1, and
  check `signal_relevance` in Stage 2.
- The 2026-06-25 example used TLEs propagated about 103 days backwards.
  Window times may be off by up to about a minute, which the padding absorbs.
  Stage-6 comparisons need contemporaneous TLEs; with the current one,
  `tart-validate` exits 3 (`TLE_STALE`, not precision validation).
- Elevations are geometric (no refraction). Skyfield's ITRS neglects polar motion (~10 m).
- Stage-1 excursions shorter than `--step` cannot be detected. Windows at the
  search-interval edge are flagged as clipped.
- Gain and phase solutions are applied as recorded and are not verified.
  Stage 3 estimates directions on each station's own timestamps, in bins on a
  common UTC grid. Stage 4 matches the nearest rows within
  `--time-match-tolerance` and never interpolates.
- Live TART APIs were not needed. At the time of writing, only the
  `za-rhodes` and `ghana` APIs answered.

## Repository layout

```text
tart-gnss/
├── README.md
├── pyproject.toml / requirements.txt
├── common_gnss.py                   # convenience wrapper for the Stage-1 CLI
├── sat_localisation.md              # Stage-1 specification
├── tart_gnss_stage2_implementation.md   # Stage-2 specification
├── tart_gnss_station_registry_implementation.md  # station-registry specification
├── tart_gnss_stage3_implementation.md   # Stage-3 specification
├── tart_gnss_stage4_implementation.md   # Stage-4 specification
├── tart_gnss_stage5_implementation.md   # Stage-5 specification
├── tart_gnss_stage6_implementation.md   # Stage-6 specification
├── tart_gnss_stage7_implementation.md   # Stage-7 specification
├── tart_gnss_campaign_analysis_implementation.md  # campaign-layer specification
├── tart_gnss_archive_scan_implementation.md       # archive-scan specification
├── tart_gnss_network_archive_discovery_implementation.md  # --all-stations specification
├── config/*.yaml                    # -> src/tart_gnss/config/*.yaml
├── configs/tart_tabsim/             # tab-sim inputs (per-station YAML, ITRF antenna files) and run records (resolved config, OMM used, logs)
├── data/stations/                   # archived TART registry snapshots + sidecars + diffs
├── data/tle/                        # archived TLE catalogues + provenance sidecars
├── data/observations/               # Stage-2 output (git-ignored)
├── data/archive_inventory/          # tart-archive-scan listing cache (git-ignored)
├── results/                         # Stage-1 example runs (2026-10-01, 2026-06-25)
├── docs/report/                     # LaTeX technical report (registry, Stages 1-2), tab-sim validation report, PDFs, figures
├── docs/tabsim_gnss/README.md       # tab-sim simulation validation: environment, commands, Ghana results
├── docs/examples/stage2/            # example Stage-2 manifests (valid and invalid)
├── docs/examples/stage3/            # synthetic end-to-end Stage-3 output
├── docs/examples/stage4/            # synthetic end-to-end Stage-4 output (with its Stage-3 input)
├── docs/examples/stage4_real/       # real 2026-06-25 Stage-4 runs (default + diagnostic) and 10 s Stage-3 input
├── docs/stage3_direction_estimation.md  # Stage-3 method, conventions, outputs
├── docs/examples/stage5/            # synthetic Stage-5 demo (Stage-4-like input, two-body/CV runs, end-to-end chain, report)
├── docs/examples/stage5_real/       # real 2026-06-25 Stage-5 runs (spec command, two-station diagnostic, CV)
├── docs/stage4_triangulation.md     # Stage-4 method, uncertainty validation, outputs
├── docs/examples/stage6/            # synthetic Stage-6 demo (known injections, ensemble, end-to-end via tart-track, report)
├── docs/examples/stage6_real/       # real 2026-06-25 Stage-6 runs (stale TLE: not precision validation)
├── docs/stage5_tracking.md          # Stage-5 method, validation, real-data results
├── docs/examples/stage7/            # synthetic Stage-7 scenarios (full Stage 3-7 chains) and demo report
├── docs/examples/stage7_real/       # real 2026-06-25 Stage-7 runs (two-station diagnostic, default)
├── docs/stage6_validation.md        # Stage-6 frames, RTN, residuals, statistics, validation, real-data results
├── docs/examples/campaign/          # synthetic 17-run campaign (runs + campaign outputs + demo report)
├── docs/examples/campaign_real/     # campaign over the real Stage-7 runs (descriptive only, exit 3)
├── docs/stage7_diagnosis.md         # Stage-7 provenance, influence, correlations, evidence rules, validation, real data
├── docs/campaign_analysis.md        # campaign inputs, cohorts, statistics, outputs, synthetic validation, real run
├── docs/examples/archive_scan_real/ # real 2026-10-09..11-03 scans (inventory-only, Stage-1 matching, two-station diagnostic)
├── docs/examples/network_archive_scan_real/ # real 2026-07-01..11-03 --all-stations scan (44 triads)
├── docs/archive_scan.md             # archive scan: method, cache, ranking, real result, decision log, handoff
├── scripts/make_report_figures.py
├── scripts/stage3_synthetic_demo.py
├── scripts/stage4_synthetic_demo.py
├── scripts/stage5_synthetic_demo.py
├── scripts/stage6_synthetic_demo.py
├── scripts/stage7_synthetic_demo.py
├── scripts/campaign_synthetic_demo.py
├── scripts/tabsim_gnss_fit.py  tabsim_gnss_diagnostics.py      # tab-sim: orbit-blind fit, MS/phase checks
├── scripts/tabsim_ghana_uncertainty.py  tabsim_ghana_rescore_truth.py  # Ghana Monte Carlo, re-scoring
├── scripts/tabsim_ghana_frame_time_audit.py  tabsim_ghana_time_scale_check.py  # frame and UTC/UT1 audits
├── src/tart_gnss/
│   ├── cli.py  stations.py  tle.py  propagation.py  visibility.py
│   ├── windows.py  search.py  output.py  provenance.py      # Stage 1
│   ├── registry.py  registry_cli.py   # TART registry: fetch, snapshot, normalise, diff
│   ├── config/  stations_registry.yaml (generated)  station_aliases.yaml  stations.yaml (legacy)
│   └── stage2/
│       ├── selection.py        # load/select Stage-1 window, query interval, signal relevance
│       ├── archive.py          # TART archive abstraction (find_observations)
│       ├── download.py         # retrieval, integrity, provenance sidecars
│       ├── native.py           # HDF5 inspection/validation, coverage
│       ├── synchronise.py      # common interval, coverage, cross-station timing
│       ├── measurement_set.py  # tart2ms wrapper + MS validation
│       ├── manifest.py         # manifest schema and writing
│       ├── errors.py           # explicit Stage-2 exceptions
│       └── cli.py              # tart-retrieve
│   └── stage3/
│       ├── loader.py           # manifest validation, MS reader, UVW check, time bins, sample selection
│       ├── geometry.py         # ENU / (l, m) / Az-El / ECEF conversions
│       ├── visibility_model.py # measurement equation, profiled chi2, Jacobian
│       ├── estimator.py        # grid search + robust least_squares
│       ├── uncertainty.py      # covariance, Az/El propagation, bootstrap
│       ├── imaging.py          # dirty-image estimator (validation)
│       ├── quality.py          # quality flags
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-localise
│   └── stage4/
│       ├── loader.py           # Stage-3 manifest/table, frozen coordinates, quality policy, sight lines
│       ├── timing.py           # time matching, angular rates, motion during offsets
│       ├── geometry.py         # Az/El and (l, m) -> ECEF lines and covariances, tangent bases
│       ├── triangulation.py    # unweighted and iteratively weighted line intersection, covariance
│       ├── uncertainty.py      # bootstrap, Monte Carlo covariance check
│       ├── quality.py          # Stage-3 policy, Stage-4 flags
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-triangulate
│   └── stage5/
│       ├── loader.py           # Stage-4 manifest/table, quality policy, covariance checks, GCRS measurements
│       ├── frames.py           # ITRS/ECEF <-> GCRS rotation of positions and covariances (Skyfield)
│       ├── dynamics.py         # constant velocity, two-body, J2; STM via variational equations; process noise
│       ├── initialisation.py   # TLE-independent initial state (GLS fit, inflated prior)
│       ├── filter.py           # EKF, Joseph update, gating, divergence recovery
│       ├── smoother.py         # RTS smoother
│       ├── innovations.py      # NIS, whiteness, run-level consistency flags
│       ├── quality.py          # Stage-4 policy, covariance integrity, Stage-5 flags
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-track
│   └── stage6/
│       ├── loader.py           # Stage-5 manifest/trajectory, full 6x6 covariance
│       ├── tle_reference.py    # --tle-file / stage1 / download selection, ensemble, provenance
│       ├── propagation.py      # SGP4 at the exact Stage-5 epochs, TLE age
│       ├── frames.py           # TEME -> GCRS chain (Skyfield), frame documentation
│       ├── rtn.py              # RTN (= RIC) basis, covariance rotation
│       ├── residuals.py        # observed - TLE residuals, Mahalanobis distances
│       ├── reference_uncertainty.py  # empirical TLE reference dispersion
│       ├── statistics.py       # summaries, trends, autocorrelation
│       ├── quality.py          # Stage-5 policy, Stage-6 flags, precision validation
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-validate
│   └── stage7/
│       ├── loader.py           # provenance chain Stage 6 -> 2 with hash verification, Stage-2 calibration lookup
│       ├── alignment.py        # time/station joins, Stage-4 minus TLE per epoch
│       ├── influence.py        # leave-one-out and pairwise re-triangulation (Stage-4 solver), direction residuals
│       ├── statistics.py       # time-detrended rank correlation, effective sample size, Holm, runs
│       ├── correlations.py     # calibration / localisation / timing / geometry / trajectory tests
│       ├── attribution.py      # station test, evidence tags and levels
│       ├── scoring.py          # LOW / MODERATE / HIGH confidence
│       ├── quality.py          # input policy, data-quality flags
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-diagnose
│   └── campaign/
│       ├── discovery.py  loader.py  schema.py  filters.py   # inputs, provenance, cohorts, exclusions
│       ├── aggregation.py       # run / epoch / station tables from completed products
│       ├── accuracy.py  station_metrics.py  constellation_metrics.py  anomaly_metrics.py  reproducibility.py
│       ├── statistics.py  bootstrap.py   # run-level statistics, comparisons, cluster bootstrap
│       ├── output.py  plotting.py  synthetic.py  errors.py
│       └── cli.py              # tart-campaign
│   └── archive_scan/
│       ├── inventory.py  cache.py       # station/day listings (Stage-2 backend), versioned checksummed cache
│       ├── intervals.py  matching.py    # filename candidate segments, intersections, Stage-1 matching
│       ├── ranking.py  output.py  errors.py
│       ├── network.py  network_run.py   # --all-stations: sites, geodesic baselines, sweep, triads
│       └── cli.py              # tart-archive-scan
└── tests/                           # 696 tests: Stage 1: 64, Stage 2: 67, registry: 36, Stage 3: 73, Stage 4: 96, Stage 5: 90, Stage 6: 62, Stage 7: 64, campaign: 58, archive scan: 83, CLI names: 3
```

## Development and testing

```bash
pip install -e ".[dev,ms,stage3]"
pytest                       # 696 tests, ~3 min, no network needed
```

The Stage-2 tests use a mocked archive and synthetic HDF5 files that follow
the real archive layout. One test runs the real `tart2ms` and is skipped if
it is not installed. The Stage-3 tests use synthetic visibilities and synthetic
Measurement Sets (written with python-casacore). `tests/stage3/test_real_data.py`
checks the visibility sign convention against the real 2026-06-25 Measurement
Sets and is skipped unless they are present (`TART_GNSS_STAGE2_DIR`).
The Stage-4 tests use synthetic Stage-3 tables built from a known 3-D truth
(exact recovery, Monte Carlo covariance validation, poor geometry, timing
offsets, a biased station, two-station and negative-range cases) and need
only numpy and pandas. The Stage-5 tests build Stage-4 tables from analytic
Kepler orbits. They check the frame rotation against Skyfield, the integrator
and STM against the analytic solution, and each of the specification's §32
scenarios. They also cover Monte Carlo NEES and coverage, the filter
lock-out, and a no-TLE run in which every TLE loader raises. The Stage-6
tests compare the TEME → GCRS chain with an independent route, and check the
RTN algebra and each of the specification's §34 cases: zero residual,
Cartesian, RTN and velocity injections, covariance rotation, TLE age, single
TLE and ensemble. They also cover TLE selection and provenance (the download
uses an injected fetcher, so no network is needed) and that the Stage-5 run
is left unmodified. The Stage-7 tests build complete Stage 3–6 chains with
injected effects: a clean network, a biased station, a calibration-gated
bias, a timing offset, poor geometry, unmodelled dynamics, a stale TLE, a
network-wide offset and missing metadata. They check the evidence tags and
levels, the leave-one-out error against truth, the provenance and hash
handling, and false-positive control on random-walk metadata. The campaign
tests build a 17-run synthetic campaign once, then check:
- that each injected population property is recovered;
- run-level bootstrap coverage and clustering;
- cohorts, filters and exclusions;
- that no pipeline stage is run during a campaign.

The archive-scan tests use a mocked archive that records every listing and
refuses downloads, together with synthetic Stage-1 JSON. They check:
- one listing per station/day;
- the edge days, and the `EMPTY`/`ERROR`/`UNKNOWN` distinctions;
- cache refresh, reuse, offline mode and checksums;
- that gaps are never filled;
- the two-station, minimum-duration and TLE-age rules;
- that ranking is not dominated by long GEO windows;
- the retrieval handoff;
- that `tart-retrieve` is never run;
- for `--all-stations`, a synthetic five-station network (one co-located
  pair, an outage, baseline-filtered triads) against hand-computed sweep
  and triad durations, a fourth station joining without fragmenting a
  triad, geodesic reference distances, and listing once per station/day.

## Technical report

[`docs/report/tart_gnss_report.pdf`](docs/report/tart_gnss_report.pdf)
(source: `tart_gnss_report.tex`).
Part I covers the TART station registry (design, snapshots, validation,
identifier resolution, change detection, the current network and its tests),
Part II covers Stage 1, Part III covers Stage 2 (archive, method,
verification and real-data results), and Part IV covers interfaces,
limitations and reproducibility. To rebuild it:

```bash
python scripts/make_report_figures.py results/common_gnss_20261001.json docs/report/figures \
    data/observations/20260625/NORAD_40534/manifest.json     # needs a --keep-native Stage-2 run
cd docs/report && latexmk -pdf tart_gnss_report.tex
```

A second, separate report,
[`docs/report/tabsim_gnss_validation.pdf`](docs/report/tabsim_gnss_validation.pdf)
(source `tabsim_gnss_validation.tex`), covers the tab-sim simulation
validation: provenance, visibility conventions, the orbit-blind fit, the
Ghana Monte Carlo and the UTC/UT1 correction. Rebuild it with
`cd docs/report && latexmk -pdf tabsim_gnss_validation.tex`.

## Authors and licence

Nadeem Oozeer, Rhodes University.

Built on [Skyfield](https://rhodesmill.org/skyfield/),
[sgp4](https://pypi.org/project/sgp4/),
[CelesTrak](https://celestrak.org) element sets, and the
[TART project](https://github.com/tart-telescope) software
([tart2ms](https://github.com/tart-telescope/tart2ms),
[tart_tools](https://pypi.org/project/tart-tools/)) and visibility archive.

No licence has been chosen yet. All rights are reserved until one is added.
