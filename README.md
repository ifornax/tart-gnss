# tart-gnss

**A pipeline for localising GNSS satellites with the TART radio-telescope network.
Stage 1 finds satellites that three TART stations can see at the same time.
Stage 2 retrieves, validates and synchronises the real TART observations for
those windows.**

| Stage | Command | What it does |
|---|---|---|
| 1 | `common-gnss` | Propagates GNSS orbits from TLEs (SGP4 via Skyfield) and finds every continuous window in which a satellite is above 20° at three stations at once. Windows are refined to 0.1 s, ranked, and stored with ITRS positions and TLE provenance. |
| 2 | `tart-retrieve` | Takes one Stage-1 window, queries the TART visibility archive, downloads the native HDF5 files with SHA-256 provenance, validates timestamps, frequency, station identity and calibration, measures the real three-station overlap, converts each station to its own Measurement Set with `tart2ms`, and writes `manifest.json` for Stage 3. |

> **Status:** Stages 1 and 2 complete, with dynamic TART station-registry integration (v0.3.0). The repository is private and under active development.

---

## Contents

- [Project context](#project-context)
- [Installation](#installation)
- [Stations](#stations)
- [Stage 1: common visibility (`common-gnss`)](#stage-1-common-visibility-common-gnss)
- [Stage 2: observation retrieval (`tart-retrieve`)](#stage-2-observation-retrieval-tart-retrieve)
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
| **1** | **Find satellites simultaneously visible from three TART stations** | **✅ `common-gnss`** |
| **2** | **Retrieve, validate and synchronise the TART data; one Measurement Set per station** | **✅ `tart-retrieve`** |
| 3 | Image the satellite independently at each station (`image-tart --manifest …`) | planned |
| 4 | Measure the satellite direction from each image | planned |
| 5 | Bayesian inference of the satellite position (x, y, z) | planned |
| 6 | Compare with the TLE-predicted position | planned |

Three conditions are kept strictly separate: **geometric visibility**
(Stage 1), **data availability** (Stage 2) and **radio detection** (Stage 3).
None is inferred from another.

Specifications: [`sat_localisation.md`](sat_localisation.md) (Stage 1) and
[`tart_gnss_stage2_implementation.md`](tart_gnss_stage2_implementation.md) (Stage 2).

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

**Reproducibility.** The live registry is for discovery only. `common-gnss`
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

## Stage 1: common visibility (`common-gnss`)

```bash
common-gnss \
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
| Stage 2, 2026-06-25, rank 8 GPS BIIF-9, 10-min subsection ([example](docs/examples/stage2/valid_20260625_NORAD_40534)) | **Valid** (exit 0). Three validated 24-antenna MSs (182 160 rows each). Real overlap 09:13:39–09:24:24 (645 s, 100 % of the retrieved interval). Station dumps are offset by 0.21–0.46 s. **9 of 24 Rhodes antennas had zero gains** (62 % flagged), and Mauritius and Madagascar each changed gain solution once. |

The example manifests were produced before the Stage-2 command was renamed,
so their recorded `command_line` starts with the old name `retrieve-tart`.
They are kept unchanged as provenance records.

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
  No epoch-matching rule between stations has been fixed yet; that is left to Stage 3.
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
├── config/*.yaml                    # -> src/tart_gnss/config/*.yaml
├── data/stations/                   # archived TART registry snapshots + sidecars + diffs
├── data/tle/                        # archived TLE catalogues + provenance sidecars
├── data/observations/               # Stage-2 output (git-ignored)
├── results/                         # Stage-1 example runs (2026-10-01, 2026-06-25)
├── docs/report/                     # LaTeX technical report (registry, Stages 1-2), PDF, figures
├── docs/examples/stage2/            # example Stage-2 manifests (valid and invalid)
├── scripts/make_report_figures.py
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
└── tests/                           # 167 tests (Stage 1: 64, Stage 2: 67, registry: 36)
```

## Development and testing

```bash
pip install -e ".[dev,ms]"
pytest                       # 167 tests, ~6 s, no network needed
```

The Stage-2 tests use a mocked archive and synthetic HDF5 files that follow
the real archive layout. One test runs the real `tart2ms` and is skipped if
it is not installed.

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
