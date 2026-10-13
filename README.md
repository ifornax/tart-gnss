# tart-gnss

**A pipeline for localising GNSS satellites with the TART radio-telescope network.
Stage 1 finds satellites that three TART stations can see at the same time.
Stage 2 retrieves, validates and synchronises the real TART observations for
those windows. Stage 3 measures the satellite's direction at each station
from the complex visibilities. Stage 4 intersects those measured lines of
sight to give the satellite's 3-D position with its covariance. Stage 5
turns those positions into a continuous position/velocity trajectory with a
full covariance, without using any TLE.**

| Stage | Command | What it does |
|---|---|---|
| 1 | `tart-common-gnss` | Propagates GNSS orbits from TLEs (SGP4 via Skyfield) and finds every continuous window in which a satellite is above 20° at three stations at once. Windows are refined to 0.1 s, ranked, and stored with ITRS positions and TLE provenance. |
| 2 | `tart-retrieve` | Takes one Stage-1 window, queries the TART visibility archive, downloads the native HDF5 files with SHA-256 provenance, validates timestamps, frequency, station identity and calibration, measures the real three-station overlap, converts each station to its own Measurement Set with `tart2ms`, and writes `manifest.json` for Stage 3. |
| 3 | `tart-localise` | Reads the Stage-2 manifest and, for each station independently, fits a point-source model directly to the complex visibilities in time bins (profiled-χ² grid search + robust `least_squares`), giving Az/El with covariance, SNR, goodness-of-fit and quality flags, cross-checked by a dirty-image peak finder. Writes CSV/Parquet tracks, `stage3_manifest.json` and diagnostic plots for Stage 4. |
| 4 | `tart-triangulate` | Reads the Stage-3 directions and the frozen station coordinates, matches stations in time (nearest within tolerance, no interpolation), and intersects the ECEF sight lines — covariance-weighted (range-scaled Stage-3 angular errors) or closed-form unweighted — giving x, y, z (ITRS/ECEF) with a 3×3 covariance, station ranges, line residuals, geometry and timing diagnostics, quality flags and an optional bootstrap. TLEs never enter the fit. |
| 5 | `tart-track` | Reads the Stage-4 positions and full 3×3 covariances, rotates them ITRS/ECEF → GCRS (Skyfield), initialises a position/velocity state from the first Stage-4 measurements only, and runs an extended Kalman filter (two-body by default; constant-velocity and J2 available) with χ²₃ innovation gating, followed by an RTS smoother. Writes the filtered and smoothed GCRS state, the ECEF position, the full 6×6 covariances, innovations, NIS/whiteness diagnostics, quality flags, `stage5_manifest.json` and seven plots. TLEs are never used. |

> **Status:** Stages 1–5 complete, with dynamic TART station-registry integration (v0.5.0). The repository is private and under active development.

---

## Contents

- [Project context](#project-context)
- [Installation](#installation)
- [Stations](#stations)
- [Stage 1: common visibility (`tart-common-gnss`)](#stage-1-common-visibility-tart-common-gnss)
- [Stage 2: observation retrieval (`tart-retrieve`)](#stage-2-observation-retrieval-tart-retrieve)
- [Stage 3: direction estimation (`tart-localise`)](#stage-3-direction-estimation-tart-localise)
- [Stage 4: 3-D localisation (`tart-triangulate`)](#stage-4-3-d-localisation-tart-triangulate)
- [Stage 5: trajectory tracking (`tart-track`)](#stage-5-trajectory-tracking-tart-track)
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
| 6 | Compare with the TLE-predicted position | planned |

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
[`docs/stage5_tracking.md`](docs/stage5_tracking.md)).

**Command names.** Every public command starts with `tart-` (`tart-stations-update`,
`tart-stations-list`, `tart-common-gnss`, `tart-retrieve`, `tart-localise`,
`tart-triangulate`, `tart-track`). This is a stable public API, and future stages must follow it;
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
within about 2%, and the coverage matches χ²₃. Method, flags and validation
are in [`docs/stage5_tracking.md`](docs/stage5_tracking.md).

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
| Stage 5, 2026-06-25 PRN 26, real data ([details](docs/stage5_tracking.md#14-real-data-run-2026-06-25-gps-prn-26-norad-40534), [example](docs/examples/stage5_real)) | **The spec command exits 3** (`INSUFFICIENT_MEASUREMENTS`), because the default Stage-4 run has one valid epoch. On the two-station Stage-4 run (35 epochs over 340 s), mean NIS is 2.8, so the *size* of the Stage-4 covariance is right. But the innovations are strongly autocorrelated (lag-1 0.69, `DYNAMICS_MISMATCH`), which shows slowly varying Stage-3 systematics, and the smoothed σ is therefore optimistic. The only three-station epoch is gated (NIS 36). Two-body and constant velocity cannot be told apart at this noise level. No TLE comparison is made (that is Stage 6). |
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
  Stage-6 comparisons need contemporaneous TLEs.
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
├── config/*.yaml                    # -> src/tart_gnss/config/*.yaml
├── data/stations/                   # archived TART registry snapshots + sidecars + diffs
├── data/tle/                        # archived TLE catalogues + provenance sidecars
├── data/observations/               # Stage-2 output (git-ignored)
├── results/                         # Stage-1 example runs (2026-10-01, 2026-06-25)
├── docs/report/                     # LaTeX technical report (registry, Stages 1-2), PDF, figures
├── docs/examples/stage2/            # example Stage-2 manifests (valid and invalid)
├── docs/examples/stage3/            # synthetic end-to-end Stage-3 output
├── docs/examples/stage4/            # synthetic end-to-end Stage-4 output (with its Stage-3 input)
├── docs/examples/stage4_real/       # real 2026-06-25 Stage-4 runs (default + diagnostic) and 10 s Stage-3 input
├── docs/stage3_direction_estimation.md  # Stage-3 method, conventions, outputs
├── docs/examples/stage5/            # synthetic Stage-5 demo (Stage-4-like input, two-body/CV runs, end-to-end chain, report)
├── docs/examples/stage5_real/       # real 2026-06-25 Stage-5 runs (spec command, two-station diagnostic, CV)
├── docs/stage4_triangulation.md     # Stage-4 method, uncertainty validation, outputs
├── docs/stage5_tracking.md          # Stage-5 method, validation, real-data results
├── scripts/make_report_figures.py
├── scripts/stage3_synthetic_demo.py
├── scripts/stage4_synthetic_demo.py
├── scripts/stage5_synthetic_demo.py
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
└── tests/                           # 428 tests: Stage 1: 64, Stage 2: 67, registry: 36, Stage 3: 73, Stage 4: 96, Stage 5: 89, CLI names: 3
```

## Development and testing

```bash
pip install -e ".[dev,ms,stage3]"
pytest                       # 428 tests, ~60 s, no network needed
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
lock-out, and a no-TLE run in which every TLE loader raises.

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

## Authors and licence

Nadeem Oozeer, Rhodes University.

Built on [Skyfield](https://rhodesmill.org/skyfield/),
[sgp4](https://pypi.org/project/sgp4/),
[CelesTrak](https://celestrak.org) element sets, and the
[TART project](https://github.com/tart-telescope) software
([tart2ms](https://github.com/tart-telescope/tart2ms),
[tart_tools](https://pypi.org/project/tart-tools/)) and visibility archive.

No licence has been chosen yet. All rights are reserved until one is added.
