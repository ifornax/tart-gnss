# Stage 1: Three-TART Common GNSS Visibility Finder

## 1. Objective

Develop a robust Python command-line tool that identifies GNSS satellites simultaneously visible from **any three selected TART radio telescope stations** over a user-specified UTC time interval.

The immediate use case is:

- TART Mauritius
- TART Madagascar
- TART South Africa
- Minimum elevation: 20 degrees
- Frequency of interest: approximately 1.575 GHz
- Satellite population: GNSS satellites

The code MUST NOT be hard-coded specifically for these three stations. The architecture must allow any three TART stations to be selected in future.

This is Stage 1 of a larger project whose eventual goal is to:

1. identify satellites simultaneously observable by three TART stations;
2. retrieve the corresponding TART visibility observations;
3. image the satellite independently at each station;
4. estimate the satellite direction from each station;
5. use Bayesian inference to infer the satellite position \((x,y,z)\);
6. compare the inferred position against the orbital position predicted from the TLE.

Stage 1 performs ONLY satellite prediction and common-visibility-window identification.

---

# 2. Scientific Definition

For satellite \(s\), station \(i\), and UTC time \(t\), calculate the topocentric elevation

\[
E_{s,i}(t).
\]

A satellite is considered visible from station \(i\) when

\[
E_{s,i}(t) \geq E_{\min},
\]

where the default threshold is

\[
E_{\min}=20^\circ.
\]

For three stations \(i=1,2,3\), define simultaneous visibility as

\[
V_s(t)=
\bigwedge_{i=1}^{3}
\left[E_{s,i}(t)\geq E_{\min}\right].
\]

A **common visibility window** is a maximal continuous interval

\[
[t_{\rm start},t_{\rm end}]
\]

during which

\[
V_s(t)=\mathrm{True}.
\]

For every such interval calculate

\[
T_{\rm common}
=
t_{\rm end}-t_{\rm start}.
\]

The primary Stage-1 result is a ranked list of satellite/common-window combinations sorted by decreasing \(T_{\rm common}\).

---

# 3. CLI Interface

Implement a CLI such as:

```bash
python common_gnss.py \
    --stations mauritius madagascar south_africa \
    --start 2026-10-01T00:00:00 \
    --end 2026-10-02T00:00:00
```

Optional arguments:

```bash
--min-elevation 20
--step 10
--tle-file gnss.tle
--output results.csv
--station-file stations.yaml
--min-duration 60
--verbose
```

Definitions:

- `--stations`
  - exactly three station identifiers for Stage 1;
  - station names must be resolved through the station configuration.

- `--start`
  - UTC ISO-8601 start time.

- `--end`
  - UTC ISO-8601 end time.

- `--min-elevation`
  - minimum topocentric elevation in degrees;
  - default = 20.

- `--step`
  - coarse orbital sampling interval in seconds;
  - default = 10 s.

- `--tle-file`
  - optional local TLE file;
  - strongly recommended for reproducibility.

- `--output`
  - output CSV filename.

- `--station-file`
  - YAML/JSON station configuration.

- `--min-duration`
  - optional minimum simultaneous visibility duration in seconds.

All internal times MUST be UTC.

Reject ambiguous timezone-naive inputs unless the documented CLI convention explicitly treats them as UTC.

---

# 4. Station Configuration

Do NOT embed station coordinates throughout the code.

Use a configuration file such as:

```yaml
stations:

  mauritius:
    name: "TART Mauritius"
    latitude_deg: null
    longitude_deg: null
    altitude_m: null

  madagascar:
    name: "TART Madagascar"
    latitude_deg: null
    longitude_deg: null
    altitude_m: null

  south_africa:
    name: "TART South Africa"
    latitude_deg: null
    longitude_deg: null
    altitude_m: null
```

IMPORTANT:

Do NOT invent coordinates.

Obtain coordinates from an authoritative TART source/API if available.

Record the source of the coordinates.

Ideally the configuration should eventually contain:

```yaml
mauritius:
  name: "TART Mauritius"
  latitude_deg: ...
  longitude_deg: ...
  altitude_m: ...
  source: "..."
```

The program should be able to support additional TART stations merely by adding entries to this file.

---

# 5. Satellite Orbit Propagation

Use TLE-based SGP4 orbit propagation.

Preferred implementation:

- `skyfield`
- `sgp4`

Do not implement orbital propagation manually.

A satellite TLE consists of:

```text
SATELLITE NAME
1 .....
2 .....
```

Instantiate each satellite using the TLE and propagate it at the requested UTC epochs.

---

# 6. GNSS Satellite Selection

The initial satellite catalogue should contain GNSS spacecraft relevant to the TART observing frequency around 1.575 GHz.

At minimum consider:

- GPS
- Galileo
- BeiDou
- GLONASS

However, do NOT assume that every spacecraft in these constellations necessarily transmits a signal relevant to 1.575 GHz.

Maintain separately:

```text
constellation
satellite_name
NORAD_ID
TLE
```

The software should support constellation filtering, for example:

```bash
--constellations GPS Galileo
```

Future extensions may distinguish actual transmitted bands/signals.

For Stage 1 the purpose is primarily geometric common visibility.

---

# 7. TLE Provenance

Scientific reproducibility is essential.

For every run record:

- TLE source
- retrieval timestamp
- TLE epoch
- satellite name
- NORAD catalogue ID
- TLE line 1
- TLE line 2

If a remote catalogue is downloaded, save the downloaded TLE catalogue locally.

For example:

```text
data/tle/
    gnss_2026-10-01T000000Z.tle
```

Never silently replace a historical TLE with a newer TLE.

If the user supplies

```bash
--tle-file ...
```

that exact catalogue must be used.

---

# 8. Important Historical-TLE Warning

Do NOT use the newest available TLE to reconstruct arbitrarily old observations without warning.

TLE/SGP4 propagation accuracy degrades as propagation moves away from the TLE epoch.

Calculate

\[
\Delta t_{\rm TLE}
=
|t_{\rm observation}-t_{\rm TLE\ epoch}|.
\]

Store this quantity or an equivalent diagnostic.

If the TLE epoch is substantially separated from the observation interval, issue a warning.

Do not invent a universal accuracy threshold unless supported by an authoritative reference.

---

# 9. Topocentric Geometry

For each station and satellite calculate:

- azimuth;
- elevation;
- range;
- optionally range rate.

At minimum retain

\[
A_i(t), E_i(t).
\]

For each satellite:

```python
station1_visible = elevation1 >= min_elevation
station2_visible = elevation2 >= min_elevation
station3_visible = elevation3 >= min_elevation

common_visible = (
    station1_visible
    & station2_visible
    & station3_visible
)
```

---

# 10. Time Sampling

Use a configurable coarse sampling interval.

Default:

```text
10 seconds
```

Do NOT assume that the coarse sample boundaries provide sufficiently accurate rise/set/common-window boundaries.

Use coarse sampling to locate transitions, then refine each boundary.

For example:

```text
19:02:20   FALSE
19:02:30   TRUE
```

implies that the threshold crossing lies between those epochs.

Refine the transition using a numerical method such as:

- bisection; or
- sufficiently fine local resampling.

Target boundary precision should be approximately 1 second unless there is a strong reason otherwise.

---

# 11. Common Visibility Windows

For each satellite identify all contiguous regions satisfying

\[
E_1\ge20^\circ,
\quad
E_2\ge20^\circ,
\quad
E_3\ge20^\circ.
\]

Do NOT simply calculate total visibility.

Preserve individual continuous windows.

Example:

```text
GPS BIIF-5

window 1:
18:22:14 – 18:47:37

window 2:
23:10:05 – 23:18:41
```

These must remain separate records.

---

# 12. Output CSV

Produce one row per:

```text
satellite × common visibility window
```

Recommended columns:

```text
rank
satellite_name
norad_id
constellation

common_start_utc
common_end_utc
duration_seconds
duration_minutes

min_elevation_deg

station_1
station_1_az_start_deg
station_1_el_start_deg
station_1_az_mid_deg
station_1_el_mid_deg
station_1_az_end_deg
station_1_el_end_deg

station_2
station_2_az_start_deg
station_2_el_start_deg
station_2_az_mid_deg
station_2_el_mid_deg
station_2_az_end_deg
station_2_el_end_deg

station_3
station_3_az_start_deg
station_3_el_start_deg
station_3_az_mid_deg
station_3_el_mid_deg
station_3_az_end_deg
station_3_el_end_deg

tle_epoch
tle_age_start_days
tle_age_mid_days
tle_source
```

Do not put multiline TLE strings directly into the primary CSV unless properly escaped.

Prefer storing complete TLE information separately.

---

# 13. Machine-Readable Detailed Output

In addition to CSV, produce JSON containing complete provenance.

Example:

```json
{
  "run": {
    "start_utc": "...",
    "end_utc": "...",
    "min_elevation_deg": 20,
    "step_seconds": 10
  },

  "stations": [
    {
      "id": "mauritius",
      "latitude_deg": 0,
      "longitude_deg": 0,
      "altitude_m": 0
    }
  ],

  "windows": [
    {
      "satellite": "...",
      "norad_id": 0,
      "constellation": "GPS",
      "start_utc": "...",
      "end_utc": "...",
      "duration_seconds": 0,
      "tle": {
        "epoch": "...",
        "line1": "...",
        "line2": "..."
      }
    }
  ]
}
```

The zeros above are schema examples ONLY.

Never substitute them for real station coordinates.

---

# 14. Position Information for Stage 2/3

The future Bayesian analysis will estimate the satellite position

\[
\mathbf r_s(t)=
\begin{bmatrix}
x(t)\\
y(t)\\
z(t)
\end{bmatrix}.
\]

Therefore Stage 1 should make it possible to recover the TLE-predicted satellite position.

For each common window midpoint, calculate or retain enough information to calculate the satellite position in an explicitly identified Earth-centred frame.

Prefer an Earth-fixed frame appropriate for comparison with terrestrial station coordinates.

Store explicitly:

```text
coordinate_frame
x_km
y_km
z_km
```

Do NOT label coordinates simply `(x,y,z)` without identifying the reference frame.

If Skyfield naturally provides another frame first, perform the appropriate transformation using established astronomy/geodesy libraries rather than implementing the transformation manually.

---

# 15. Future Bayesian Requirement

The Stage-1 design must anticipate the later model.

For station \(i\), with Earth-fixed location

\[
\mathbf r_i,
\]

and satellite location

\[
\mathbf r_s,
\]

the geometric line of sight is

\[
\mathbf d_i
=
\mathbf r_s-\mathbf r_i.
\]

The unit line-of-sight vector is

\[
\hat{\mathbf d}_i
=
\frac{\mathbf d_i}
{\|\mathbf d_i\|}.
\]

Later, TART imaging will provide measured directions

\[
(A_i^{\rm obs},E_i^{\rm obs})
\]

or equivalent celestial coordinates.

The Bayesian model will compare observed directions against those predicted from a candidate satellite position.

Therefore preserve:

- exact station coordinates;
- exact UTC epochs;
- TLE provenance;
- predicted azimuth/elevation;
- predicted satellite position.

---

# 16. Software Structure

Do NOT write everything in one script.

Recommended structure:

```text
tart_gnss/
│
├── README.md
├── pyproject.toml
├── requirements.txt
│
├── config/
│   └── stations.yaml
│
├── data/
│   └── tle/
│
├── src/
│   └── tart_gnss/
│       ├── __init__.py
│       ├── cli.py
│       ├── stations.py
│       ├── tle.py
│       ├── propagation.py
│       ├── visibility.py
│       ├── windows.py
│       ├── output.py
│       └── provenance.py
│
└── tests/
    ├── test_stations.py
    ├── test_tle.py
    ├── test_visibility.py
    ├── test_windows.py
    └── test_cli.py
```

---

# 17. Module Responsibilities

## `stations.py`

Responsible for:

- loading station configuration;
- validating latitude;
- validating longitude;
- validating altitude;
- resolving station IDs.

Reject unknown station names with a clear error.

---

## `tle.py`

Responsible for:

- reading local TLE files;
- parsing satellite names;
- obtaining NORAD IDs;
- constellation classification;
- TLE provenance.

Network download logic, if implemented, must be isolated here.

---

## `propagation.py`

Responsible for:

- SGP4 propagation;
- topocentric az/el;
- satellite range;
- satellite position.

No file-output logic should exist here.

---

## `visibility.py`

Responsible for evaluating

\[
E_i(t)\ge E_{\min}.
\]

Keep this logic generic so it can eventually support \(N\) stations.

Although Stage 1 CLI requires exactly three stations, internal code should preferably operate on:

```python
stations: Sequence[Station]
```

rather than:

```python
station1
station2
station3
```

This will make future 2-, 3-, 4-, or N-station experiments straightforward.

---

## `windows.py`

Responsible for:

- detecting FALSE→TRUE transitions;
- detecting TRUE→FALSE transitions;
- refining boundaries;
- creating continuous common-visibility windows;
- calculating durations.

---

## `output.py`

Responsible for:

- CSV;
- JSON;
- sorting;
- ranking.

---

# 18. Logging

Use Python's `logging` module.

Example:

```python
logging.info(
    "Searching %d GNSS satellites between %s and %s",
    len(satellites),
    start,
    end,
)
```

Log:

- stations selected;
- station coordinates;
- number of satellites;
- TLE source;
- time interval;
- sampling interval;
- elevation threshold;
- number of common windows;
- output filenames.

Do not use `print()` for operational diagnostics.

CLI result summaries may use normal terminal output if appropriate.

---

# 19. Docstrings

Every public function and class MUST have a clear docstring describing:

- purpose;
- inputs;
- outputs;
- units;
- coordinate system where relevant;
- important assumptions.

Example:

```python
def compute_topocentric_position(...):
    """Compute satellite azimuth and elevation for a ground station.

    Parameters
    ----------
    satellite
        Skyfield EarthSatellite object.
    station
        Ground-station location.
    times
        UTC Skyfield time array.

    Returns
    -------
    azimuth_deg : numpy.ndarray
        Topocentric azimuth in degrees.
    elevation_deg : numpy.ndarray
        Topocentric elevation in degrees.
    range_km : numpy.ndarray
        Observer-to-satellite range in kilometres.

    Notes
    -----
    Visibility is determined separately using the configured minimum
    elevation threshold.
    """
```

---

# 20. Tests

Tests are mandatory.

At minimum test:

### Station validation

Invalid:

```text
latitude = 100°
longitude = 300°
```

must fail.

### Visibility logic

Synthetic elevations:

```text
station A = [19, 21, 25, 18]
station B = [22, 23, 24, 25]
station C = [30, 30, 19, 30]
```

with threshold 20° should produce:

```text
False
True
False
False
```

### Window extraction

Given:

```text
False False True True True False
```

verify exactly one window.

### Multiple windows

Given:

```text
False True True False False True True False
```

verify exactly two independent windows.

### Boundary refinement

Verify that the refined threshold crossing falls within the original coarse bracket.

### CLI

Invalid dates such as:

```text
--start later-than-end
```

must fail clearly.

Unknown stations must fail clearly.

---

# 21. Reproducibility Metadata

Every run should produce something like:

```text
results/
    common_gnss_20261001.csv
    common_gnss_20261001.json
    run_metadata_20261001.json
```

Metadata should include:

```text
software version
Python version
Skyfield version
SGP4 version
run timestamp
start
end
station IDs
station coordinates
minimum elevation
sampling interval
TLE source
TLE filename
TLE retrieval timestamp
```

If possible also record a SHA-256 hash of the TLE input file.

---

# 22. Stage-1 Terminal Summary

At completion display a concise table:

```text
COMMON GNSS VISIBILITY

Stations:
  Mauritius
  Madagascar
  South Africa

Interval:
  2026-10-01T00:00:00Z
  2026-10-02T00:00:00Z

Minimum elevation: 20 deg

Rank Satellite       NORAD  Start UTC             End UTC               Duration
1    GPS ...         .....  2026-10-01T...        2026-10-01T...        38.4 min
2    Galileo ...     .....  2026-10-01T...        2026-10-01T...        31.2 min
...
```

Do not fabricate results if no satellite satisfies the condition.

Instead report:

```text
No common GNSS visibility windows found.
```

---

# 23. Do Not Implement Yet

Stage 1 MUST NOT yet implement:

- TART visibility download;
- `tart2ms` conversion;
- Measurement Set processing;
- WSClean imaging;
- source extraction from TART images;
- Bayesian inference;
- MCMC;
- nested sampling;
- triangulation from measured TART directions.

Those belong to later stages.

However, Stage 1 outputs MUST contain sufficient timing, station, TLE, and coordinate information to support them.

---

# 24. Stage 2 Interface

Design Stage 1 so that its output can eventually drive:

```bash
python retrieve_tart.py \
    --visibility-results common_gnss.csv \
    --rank 1
```

This will eventually:

1. read the selected common visibility window;
2. determine the three TART station APIs/archive names;
3. query their archives;
4. retrieve the corresponding HDF5 visibility data;
5. use `tart2ms` to generate three Measurement Sets.

The current `tart2ms` implementation supports HDF5 archive visibility data and conversion into Measurement Sets. Its archive-search functionality should be investigated and reused rather than independently reimplementing the entire TART archive client.

---

# 25. Acceptance Criteria

Stage 1 is complete when the following works:

```bash
python common_gnss.py \
    --stations mauritius madagascar south_africa \
    --start 2026-10-01T00:00:00Z \
    --end 2026-10-02T00:00:00Z \
    --min-elevation 20 \
    --step 10
```

and produces:

1. validated three-station configuration;
2. GNSS TLE catalogue with provenance;
3. SGP4 propagation;
4. az/el calculation for each station;
5. simultaneous elevation >20° detection;
6. refined continuous common windows;
7. ranked common-duration results;
8. CSV output;
9. detailed JSON output;
10. TLE provenance;
11. predicted satellite coordinates with explicit reference frame;
12. reproducibility metadata;
13. passing automated tests.

The implementation must prioritize scientific reproducibility and correctness over code brevity.

Do not invent missing station coordinates, TLE information, API endpoints, satellite properties, or observational data. If authoritative information cannot be obtained, fail explicitly or require it from configuration.
