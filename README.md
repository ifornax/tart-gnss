# tart-gnss

**Stage 1 of the TART GNSS satellite-localisation pipeline: finding GNSS
satellites that several TART radio telescopes can see at the same time.**

`tart_gnss` propagates GNSS orbits from TLEs (SGP4 via Skyfield). It finds every
continuous time window in which a satellite is above a minimum elevation
(default 20°) at **three TART stations at once**, refines each window's
boundaries to 0.1 s and ranks the windows by duration. Each window carries
everything later stages need: exact UTC epochs, station coordinates, predicted
azimuth/elevation, the satellite position in an Earth-fixed frame (ITRS) and
full TLE provenance.

> **Status:** Stage 1 complete (v0.1.0). Private repository, under active development.

---

## Contents

- [Project context](#project-context)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command-line options](#command-line-options)
- [Outputs](#outputs)
- [Stations](#stations)
- [TLEs and reproducibility](#tles-and-reproducibility)
- [Method in brief](#method-in-brief)
- [Caveats](#caveats)
- [Repository layout](#repository-layout)
- [Development and testing](#development-and-testing)
- [Technical report](#technical-report)
- [Authors and licence](#authors-and-licence)

---

## Project context

The [Transient Array Radio Telescope (TART)](https://tart.elec.ac.nz) is a
low-cost, open-source 24-element all-sky interferometer operating at the GPS L1
frequency (1575.42 MHz). Three TART stations that image the same satellite at
the same time give three lines of sight, so the satellite's 3-D position can
be inferred and compared with its orbital prediction.

| Stage | Task | Status |
|---|---|---|
| **1** | **Find satellites simultaneously visible from three TART stations** | **✅ this repository** |
| 2 | Retrieve the matching TART visibilities (reusing `tart2ms`) | planned |
| 3 | Image the satellite independently at each station | planned |
| 4 | Measure the satellite direction from each image | planned |
| 5 | Bayesian inference of the satellite position (x, y, z) | planned |
| 6 | Compare with the TLE-predicted position | planned |

The full Stage-1 specification is in [`sat_localisation.md`](sat_localisation.md).

## Installation

Requires Python ≥ 3.9. Dependencies (NumPy, Skyfield, sgp4, PyYAML) are installed automatically.

```bash
git clone git@github.com:ifornax/tart-gnss.git
cd tart-gnss
pip install .                 # regular install
pip install -e ".[dev]"       # editable, with pytest and matplotlib
```

Because the repository is private, a direct install needs SSH access:

```bash
pip install "git+ssh://git@github.com/ifornax/tart-gnss.git"
```

## Quick start

```bash
common-gnss \
    --stations mauritius madagascar south_africa \
    --start 2026-10-01T00:00:00Z \
    --end   2026-10-02T00:00:00Z \
    --min-elevation 20 --step 10
```

`python common_gnss.py …` and `python -m tart_gnss …` are equivalent.

Example terminal summary (abridged):

```text
COMMON GNSS VISIBILITY

Stations:
  TART Mauritius (Bel Air) [mauritius]
  TART Madagascar (Antananarivo) [madagascar]
  TART South Africa (Rhodes University, Makhanda) [south_africa]

Interval:
  2026-10-01T00:00:00Z
  2026-10-02T00:00:00Z

Minimum elevation: 20 deg

Rank Satellite                   NORAD  Const.   Start UTC             End UTC                 Duration
-------------------------------------------------------------------------------------------------------
1    BEIDOU-2 G5 (C05)           38091  BeiDou   2026-10-01T00:00:00Z  2026-10-02T00:00:00Z    1440.0 min*
3    BEIDOU-2 G6 (C02)           38953  BeiDou   2026-10-01T12:28:02Z  2026-10-01T20:16:28Z     468.4 min
4    GSAT0219 (GALILEO 23)       43566  Galileo  2026-10-01T08:28:07Z  2026-10-01T15:46:19Z     438.2 min
...
* window truncated by the search interval
Total: 154 common visibility windows
```

If no satellite qualifies, the tool prints `No common GNSS visibility windows found.`

## Command-line options

| Option | Meaning | Default |
|---|---|---|
| `--stations` | exactly three station IDs from the station file | required |
| `--start`, `--end` | ISO-8601 times. `Z`/offsets are honoured; **naive times are interpreted as UTC** | required |
| `--min-elevation` | elevation threshold, degrees | 20 |
| `--step` | coarse sampling step, seconds | 10 |
| `--tolerance` | boundary refinement tolerance, seconds | 0.1 |
| `--tle-file` | use this TLE catalogue exactly | download from CelesTrak |
| `--tle-dir` | where downloaded catalogues are archived | `data/tle` |
| `--constellations` | restrict to e.g. `GPS Galileo` (`GPS`, `Galileo`, `BeiDou`, `GLONASS`) | all |
| `--output` | CSV path | `results/common_gnss_<YYYYMMDD>.csv` |
| `--station-file` | YAML/JSON station configuration | bundled file |
| `--min-duration` | drop windows shorter than this, seconds | 0 |
| `--tle-age-warning-days` | warn when \|t_obs − t_TLE epoch\| exceeds this (a flag level, not an accuracy bound) | 3 |
| `--summary-rows` | rows shown in the terminal summary | 20 |
| `--verbose` | debug logging | off |

Input errors (unknown station, malformed date, start not before end, wrong
number of stations) exit with status 2 and a clear message.

## Outputs

Each run writes three files (shown for the default name):

| File | Contents |
|---|---|
| `results/common_gnss_20261001.csv` | One row per satellite × window: rank, satellite, NORAD ID, constellation, start/end UTC, duration, clipping flags, and az/el at start/mid/end plus mid-window range for each station. Also the mid-window satellite position in **ITRS** (km), TLE epoch, TLE age and TLE source. |
| `results/common_gnss_20261001.json` | Everything above plus full TLE lines, range rates, ITRS position and velocity at start/mid/end, per-station line-of-sight unit vectors, station ITRS vectors and the list of every searched satellite. |
| `results/run_metadata_20261001.json` | Software versions (tart_gnss, Python, Skyfield, sgp4, NumPy), command line, stations, thresholds, TLE file with SHA-256 and retrieval time, station-file SHA-256 and a TLE-age summary. |

All times are UTC (`…Z`). Cartesian coordinates are always labelled with their
frame (`coordinate_frame = "ITRS"`).

## Stations

Coordinates are never hard-coded. They live in
[`src/tart_gnss/config/stations.yaml`](src/tart_gnss/config/stations.yaml)
(symlinked as `config/stations.yaml`). They were taken from the TART registry
`https://map.elec.ac.nz/api/v1/telescopes` on 2026-10-07, and each entry
records its source.

| ID | TART id | Site | Lat (°) | Lon (°) | Alt (m) |
|---|---|---|---|---|---|
| `mauritius` | `mu-udm` | Bel Air, Mauritius | −20.2587508 | 57.7591989 | 20.0 |
| `madagascar` | `mg-tana` | Antananarivo, Madagascar | −18.8971446 | 47.5551677 | 1280.0 |
| `south_africa` | `za-rhodes` | Rhodes University, Makhanda | −33.319413 | 26.507688 | 642.26 |
| `south_africa_hammanskraal` | `za-hammanskraal` | Kwalata, Hammanskraal | −25.38948 | 28.32231 | 1163.0 |
| `namibia` | `na-unam` | UNAM, Windhoek | −22.612053 | 17.856784 | 1655.0 |
| `zambia` | `zm-cbu` | Copperbelt University | −12.80527778 | 28.23722222 | 1270.0 |
| `ghana` | `ghana` | GRAO | 5.749725 | −0.305557 | 72.0 |

To add a station, add an entry with `latitude_deg`, `longitude_deg`,
`altitude_m` and `source`; no code change is needed. Missing (`null`),
non-numeric or out-of-range values are rejected. These are the coordinates
the operators configured on each telescope. Their survey accuracy is
undocumented, so use surveyed values before Stage 5.

## TLEs and reproducibility

- Without `--tle-file`, current GNSS TLEs are downloaded from CelesTrak
  (`gps-ops`, `galileo`, `beidou`, `glo-ops`). They are archived as
  `data/tle/gnss_<retrieval time>.tle` with a `.meta.json` sidecar holding
  the source URLs, retrieval time, SHA-256 and the constellation of each NORAD ID.
- Archived catalogues are **never overwritten**. If a file no longer matches
  its recorded SHA-256, the run stops with an error.
- **CelesTrak serves only current elements.** To analyse an older observation,
  supply a catalogue from the observation date with `--tle-file`. The tool
  warns when satellites are propagated backwards from their TLE epoch, or
  further than `--tle-age-warning-days` from it.
- To reproduce the example results in this repository:

  ```bash
  common-gnss --stations mauritius madagascar south_africa \
      --start 2026-10-01T00:00:00Z --end 2026-10-02T00:00:00Z \
      --tle-file data/tle/gnss_2026-10-07T101156Z.tle
  ```

## Method in brief

1. Propagate each satellite with SGP4 (Skyfield) on a coarse grid (default 10 s).
2. Compute the geometric topocentric elevation `E_i(t)` at every station.
3. Common visibility is `V(t) = AND_i [E_i(t) ≥ E_min]`, equivalently `min_i E_i(t) ≥ E_min`.
4. Split `V` into maximal runs. Refine each rise or set by bisection on the
   exact predicate to within 0.1 s, so the refined boundary always stays inside
   the original coarse bracket.
5. For each window, record the az/el/range at each station at start, middle and
   end, the satellite's ITRS position and velocity, and the line-of-sight
   vectors `d_i = r_s − r_i`.
6. Rank windows by decreasing duration. Ties are broken by start time, then NORAD ID.

## Caveats

- **Geometric visibility ≠ detectability.** TART records 2.5 MHz around
  1575.42 MHz. GPS L1, Galileo E1 and BeiDou-3 B1C fall in this band, but
  BeiDou-2 (B1I, 1561 MHz) and legacy GLONASS FDMA (1598–1606 MHz) do not.
  For Stage 2, prefer `--constellations GPS Galileo` or BeiDou-3 windows.
- Elevations are geometric: refraction is not applied (about 2.7′ at 20°).
- Skyfield's ITRS neglects polar motion (~10 m).
- Excursions shorter than `--step` cannot be detected.
- Windows touching the search-interval edge are flagged as clipped.
- At the time of writing, only the `za-rhodes` and `ghana` TART APIs responded.
  Verify the others before Stage 2.

## Repository layout

```text
tart-gnss/
├── README.md
├── pyproject.toml / requirements.txt
├── common_gnss.py               # convenience wrapper for the CLI
├── sat_localisation.md          # Stage-1 specification
├── config/stations.yaml         # -> src/tart_gnss/config/stations.yaml
├── data/tle/                    # archived TLE catalogues + provenance sidecars
├── results/                     # example run (2026-10-01)
├── docs/report/                 # LaTeX technical report, PDF and figures
├── scripts/make_report_figures.py
├── src/tart_gnss/
│   ├── cli.py                   # argument parsing, validation, logging
│   ├── stations.py              # station configuration and validation
│   ├── tle.py                   # TLE parsing, provenance, CelesTrak download
│   ├── propagation.py           # SGP4, az/el/range, ITRS positions
│   ├── visibility.py            # N-station elevation-threshold logic
│   ├── windows.py               # runs and bisection refinement
│   ├── search.py                # orchestration and window characterisation
│   ├── output.py                # ranking, CSV, JSON, terminal summary
│   └── provenance.py            # versions, hashes, run metadata
└── tests/                       # pytest suite (62 tests)
```

## Development and testing

```bash
pip install -e ".[dev]"
pytest                       # 62 tests, < 1 s, no network needed
```

The tests cover station validation, TLE parsing and provenance, visibility
logic, window extraction and refinement, and CLI error handling. They also
check propagation physics: orbital radii, and agreement better than 0.001°
between the az/el results and an independent ENU projection of the ITRS
line of sight.

## Technical report

[`docs/report/tart_gnss_stage1.pdf`](docs/report/tart_gnss_stage1.pdf)
(source: `tart_gnss_stage1.tex`) covers the problem definition, data
provenance, method, software design, verification, the 2026-10-01 results and
the interface to Stages 2–6. To rebuild it:

```bash
python scripts/make_report_figures.py results/common_gnss_20261001.json docs/report/figures
cd docs/report && latexmk -pdf tart_gnss_stage1.tex
```

## Authors and licence

Nadeem Oozeer, Rhodes University.

Built on [Skyfield](https://rhodesmill.org/skyfield/),
[sgp4](https://pypi.org/project/sgp4/),
[CelesTrak](https://celestrak.org) element sets and the
[TART project](https://github.com/tart-telescope).

No licence has been chosen yet. All rights are reserved until one is added.
