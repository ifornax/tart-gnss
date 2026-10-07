# tart-gnss — Stage 1: common GNSS visibility for TART stations

`tart_gnss` finds the GNSS satellites (GPS, Galileo, BeiDou, GLONASS) that are
**simultaneously above a minimum elevation at three TART radio-telescope
stations**, and ranks the continuous common-visibility windows by duration.

It is Stage 1 of a pipeline whose goal is to localise satellites in 3-D from
independent TART images (retrieve visibilities → image each station →
measure directions → Bayesian position inference → compare with the TLE
prediction). Only Stage 1 is implemented here; its outputs carry the station,
timing, TLE and Earth-fixed coordinate information the later stages need.

The full specification is in [`sat_localisation.md`](sat_localisation.md) and
the technical report in [`docs/report/tart_gnss_stage1.pdf`](docs/report/tart_gnss_stage1.pdf).

## Installation

```bash
pip install .                 # from a checkout
pip install -e ".[dev]"       # editable, with pytest + matplotlib
# or directly from GitHub
pip install git+https://github.com/ifornax/tart-gnss.git
```

Requires Python ≥ 3.9, NumPy, Skyfield, sgp4 and PyYAML.

## Usage

```bash
common-gnss \
    --stations mauritius madagascar south_africa \
    --start 2026-10-01T00:00:00Z \
    --end   2026-10-02T00:00:00Z \
    --min-elevation 20 --step 10
```

`python common_gnss.py ...` and `python -m tart_gnss ...` are equivalent.

| Option | Meaning (default) |
|---|---|
| `--stations` | exactly three station IDs from the station file |
| `--start`, `--end` | ISO-8601 UTC. `Z`/offsets honoured; **naive times are interpreted as UTC** |
| `--min-elevation` | elevation threshold, deg (20) |
| `--step` | coarse sampling step, s (10) |
| `--tolerance` | boundary refinement tolerance, s (0.1) |
| `--tle-file` | use this TLE catalogue exactly (recommended for reproducibility) |
| `--tle-dir` | where downloaded catalogues are archived (`data/tle`) |
| `--constellations` | e.g. `GPS Galileo` |
| `--output` | CSV path (`results/common_gnss_<YYYYMMDD>.csv`) |
| `--station-file` | YAML/JSON station configuration (bundled file) |
| `--min-duration` | drop windows shorter than this, s (0) |
| `--tle-age-warning-days` | warn when \|t_obs − t_TLE\| exceeds this (3; a flag level, not an accuracy bound) |
| `--verbose` | debug logging |

Without `--tle-file`, current TLEs are downloaded from CelesTrak and archived
as `data/tle/gnss_<retrieval>.tle` plus a `.meta.json` sidecar (source URLs,
retrieval time, SHA-256, constellation per NORAD ID). Archived files are never
overwritten. **CelesTrak only serves current elements**: to analyse an old
observation, supply a contemporaneous catalogue with `--tle-file`; the tool
warns when satellites are propagated backwards or far from their TLE epoch.

### Outputs

* `common_gnss_<tag>.csv` — one row per satellite × window: rank, times,
  duration, clipping flags, az/el at start/mid/end for each station, mid-window
  range, mid-window satellite position in **ITRS** (km), TLE epoch/age/source.
* `common_gnss_<tag>.json` — the same plus full TLE lines, range rates,
  ITRS position/velocity at start/mid/end, line-of-sight unit vectors per
  station, station ITRS vectors, and every searched satellite.
* `run_metadata_<tag>.json` — software versions, command line, stations,
  thresholds, TLE file + SHA-256 + retrieval time, TLE-age summary.

## Stations

Station coordinates live only in
[`src/tart_gnss/config/stations.yaml`](src/tart_gnss/config/stations.yaml)
(symlinked as `config/stations.yaml`). They were taken from the TART registry
`https://map.elec.ac.nz/api/v1/telescopes` on 2026-10-07; each entry records
its source. `south_africa` is the Rhodes University telescope (`za-rhodes`);
Hammanskraal, Namibia, Zambia and Ghana are also configured. Add a station by
adding an entry — no code change is needed. Missing (`null`) or out-of-range
coordinates are rejected.

## Caveats

* Visibility is purely geometric. BeiDou-2 (B1I, 1561 MHz) and legacy
  GLONASS FDMA (1598–1606 MHz) do not transmit within TART's 2.5 MHz band at
  1575.42 MHz; prefer `--constellations GPS Galileo` (or BeiDou-3) for Stage 2.
* Elevations are geometric (no refraction); ITRS ignores polar motion (~10 m).
* Windows touching the search interval edge are flagged as clipped.

## Development

```bash
pip install -e ".[dev]"
pytest                                     # 62 tests
python scripts/make_report_figures.py results/common_gnss_20261001.json docs/report/figures
cd docs/report && latexmk -pdf tart_gnss_stage1.tex
```

Package layout: `stations.py` (configuration), `tle.py` (TLEs, provenance,
download), `propagation.py` (SGP4, az/el, ITRS), `visibility.py`,
`windows.py` (runs + bisection refinement), `search.py` (orchestration),
`output.py` (CSV/JSON/summary), `provenance.py`, `cli.py`.
