# Stage 2 --- TART Observation Retrieval, Synchronisation, and Measurement-Set Preparation

## 1. Objective

Implement Stage 2 of the `tart-gnss` pipeline.

Stage 1 identifies GNSS satellites simultaneously geometrically visible
from three selected TART stations. Stage 2 uses a selected Stage-1
common-visibility window to locate, retrieve, validate, and prepare the
corresponding **real TART interferometric observations** from all three
stations.

The processing chain is:

``` text
Stage-1 common GNSS window
          |
          v
Select satellite/window
          |
          v
Resolve TART station/archive identifiers
          |
          v
Query TART archives
          |
          v
Retrieve native visibility files
          |
          v
Validate timestamps and metadata
          |
          v
Determine actual 3-station temporal overlap
          |
          v
Convert each station to Measurement Set
          |
          v
Stage-2 observation manifest
          |
          v
Stage 3: independent imaging
```

Stage 2 MUST NOT perform satellite localisation or Bayesian inference.

## 2. Scientific Requirement

For a Stage-1 common visibility interval

$$[t_{pred,start},t_{pred,end}],$$

retrieve observations from stations $i=1,2,3$.

The usable three-station interval is

$$t_{common,start}=\max_i(t_{i,start})$$

and

$$t_{common,end}=\min_i(t_{i,end}).$$

A valid three-station dataset requires

$$t_{common,end}>t_{common,start}.$$

Do not call an observation a three-station observation unless all three
stations contain real visibility data covering a non-zero common
interval.

## 3. Stage-1 Input

Prefer the detailed Stage-1 JSON output:

``` bash
retrieve-tart \
  --visibility-results results/common_gnss_20261001.json \
  --rank 1
```

Support selection by `--rank` and `--norad-id`, with optional
`--window-index`.

Preserve: - satellite name, NORAD ID and constellation; - common
start/end UTC; - three station IDs and coordinates; - predicted az/el
and ITRS position; - TLE epoch and provenance.

Do not silently recompute Stage-1 geometry.

## 4. Satellite Selection

For initial experiments prioritize signals inside the TART band around
1575.42 MHz:

1.  GPS L1;
2.  Galileo E1;
3.  BeiDou-3 B1C where appropriate.

Geometric visibility and radio detectability are separate conditions.
Record a `signal_relevance` field such as `GPS_L1`, `GALILEO_E1`,
`BEIDOU_B1C`, or `UNKNOWN`. Do not fabricate assignments.

## 5. Retrieval Padding

Support configurable padding:

``` bash
--padding 120
```

with

$$t_{query,start}=t_{pred,start}-P$$

and

$$t_{query,end}=t_{pred,end}+P.$$

Default: 120 seconds.

## 6. CLI

Recommended interface:

``` bash
retrieve-tart \
  --visibility-results results/common_gnss_20261001.json \
  --rank 1 \
  --padding 120 \
  --output-dir data/observations
```

Arguments:

``` text
--visibility-results PATH
--rank INTEGER
--norad-id INTEGER
--window-index INTEGER
--padding SECONDS
--output-dir PATH
--download
--convert-ms
--keep-native
--allow-partial
--dry-run
--verbose
```

Selection arguments must identify one unambiguous window.

## 7. TART Station Mapping

Extend the existing station configuration rather than hard-coding
archive identifiers in retrieval functions.

Example:

``` yaml
stations:
  mauritius:
    name: "TART Mauritius"
    tart_id: "mu-udm"
    latitude_deg: ...
    longitude_deg: ...
    altitude_m: ...
    api_url: null
    archive_name: null
    source: "TART registry"
```

Determine actual API/archive identifiers from authoritative TART
sources. Do not invent URLs or archive names.

## 8. Reuse Existing TART Software

Reuse `tart2ms`, `tart_tools`, TART archive handlers, and maintained
TART API clients where practical.

Do not independently reimplement archive communication when maintained
TART software already provides it.

`tart2ms` should be the preferred path for conversion of native TART
visibility data into CASA Measurement Sets.

## 9. Archive Discovery

Implement a normalized archive-discovery abstraction:

``` python
def find_observations(station, start_utc, end_utc):
    """Find TART observations intersecting a UTC interval.

    Parameters
    ----------
    station
        Validated TART station configuration.
    start_utc, end_utc
        UTC query boundaries.

    Returns
    -------
    list
        Normalized observation records containing archive identifiers,
        filenames, timestamps and available format metadata.
    """
```

Return, where available:

``` text
station_id
tart_id
filename
remote_identifier
start_utc
end_utc
format
size_bytes
download/archive reference
```

Do not expose raw HTTP responses to the rest of the pipeline.

## 10. Missing Stations

If one station has no data:

``` text
Mauritius      FOUND
Madagascar     NOT FOUND
South Africa   FOUND
```

the program MUST NOT silently label the result a three-station
observation.

Set:

``` text
three_station_valid: false
```

`--allow-partial` may permit diagnostic processing, but partial data
must remain explicitly marked invalid for three-station localisation.

## 11. Native Data

Prefer retaining original TART archive files. HDF5 should be preferred
where available for multi-timestamp observations.

Recommended layout:

``` text
data/observations/
└── 20261001/
    └── NORAD_43566/
        ├── stage1_window.json
        ├── manifest.json
        ├── mauritius/raw/
        ├── madagascar/raw/
        └── south_africa/raw/
```

Never silently overwrite existing files.

## 12. File Integrity and Provenance

For every native file record:

-   SHA-256;
-   file size;
-   remote identifier;
-   retrieval UTC;
-   original filename;
-   local filename.

Reuse an existing local file only when its provenance/integrity is
consistent.

## 13. Native Observation Validation

Inspect the actual data, not merely filenames.

At minimum validate:

-   timestamps;
-   station identity;
-   number of antennas where available;
-   number of baselines;
-   number of time samples;
-   observing frequency;
-   visibility presence;
-   calibration metadata where available;
-   flags where available.

## 14. Frequency Validation

Record actual observing frequency and bandwidth from the observation
metadata.

Verify consistency with the expected TART band around 1575.42 MHz. Do
not invent an acceptable tolerance before inspecting the instrumental
metadata.

If frequency metadata are unavailable, report that explicitly.

## 15. Temporal Coverage

For every station calculate:

``` text
first_timestamp_utc
last_timestamp_utc
number_of_timestamps
median_time_spacing_seconds
```

Detect significant gaps. File presence alone does not establish
continuous temporal coverage.

## 16. Three-Station Synchronisation

Calculate:

``` text
actual_common_start_utc
actual_common_end_utc
actual_common_duration_seconds
```

Characterize cross-station timing rather than assuming exact dump
synchronization.

Record:

``` text
station_1 median cadence
station_2 median cadence
station_3 median cadence
median nearest timestamp separation
maximum nearest timestamp separation
```

The final epoch-matching rule for localisation can be refined after real
data are inspected.

## 17. Coverage of Predicted Window

For each station and for the three-station intersection calculate the
fraction of the Stage-1 predicted interval actually covered by data:

$$
f_{coverage} =
\frac{T(\mathrm{actual\ data}\cap\mathrm{predicted\ window})}
{T_{predicted}}.
$$

Record:

``` text
predicted_window_duration_s
station_1_coverage_fraction
station_2_coverage_fraction
station_3_coverage_fraction
three_station_coverage_fraction
```

## 18. Measurement Set Conversion

Create one independent Measurement Set per TART station:

``` text
mauritius.ms
madagascar.ms
south_africa.ms
```

Do NOT concatenate stations into a single interferometer. They are
independent arrays whose images will later supply independent lines of
sight.

## 19. tart2ms Wrapper

Use a well-defined wrapper, for example:

``` python
def convert_to_ms(input_files, output_ms, station):
    """Convert native TART visibility files into one station Measurement Set.

    Parameters
    ----------
    input_files
        Ordered native TART observation files.
    output_ms
        Destination Measurement Set.
    station
        Station metadata used for provenance and validation.

    Returns
    -------
    object
        Validated Measurement-Set metadata.

    Notes
    -----
    The function must capture conversion failures and verify the generated
    Measurement Set before returning successfully.
    """
```

Prefer a Python API or safe subprocess argument list over shell-string
construction.

## 20. Measurement Set Validation

After conversion verify:

``` text
MS exists
non-zero rows
TIME readable
UVW present
DATA/visibility column present
ANTENNA table present
SPECTRAL_WINDOW present
```

Record:

``` text
number_of_rows
number_of_antennas
number_of_channels
frequency
first_time
last_time
```

## 21. Calibration

Stage 2 must not invent new calibration algorithms.

Preserve existing calibration information and document how `tart2ms`
handles it.

Record:

``` text
calibration_status
calibration_source
```

Do not claim observations are calibrated unless verified.

## 22. Output Layout

``` text
data/observations/
└── 20261001/
    └── NORAD_43566/
        ├── stage1_window.json
        ├── manifest.json
        ├── retrieval.log
        ├── mauritius/
        │   ├── raw/
        │   └── ms/mauritius.ms
        ├── madagascar/
        │   ├── raw/
        │   └── ms/madagascar.ms
        └── south_africa/
            ├── raw/
            └── ms/south_africa.ms
```

## 23. Stage-2 Manifest

`manifest.json` is the primary Stage-2 → Stage-3 interface.

Recommended structure:

``` json
{
  "stage": 2,
  "three_station_valid": true,
  "satellite": {},
  "prediction": {},
  "stations": {},
  "synchronisation": {},
  "software": {},
  "provenance": {}
}
```

### Satellite section

Include:

``` text
satellite_name
norad_id
constellation
signal_relevance
tle_epoch
tle_source
```

### Prediction section

Include:

``` text
predicted_common_start_utc
predicted_common_end_utc
predicted_duration_seconds
query_start_utc
query_end_utc
padding_seconds
stage1_result_file
stage1_result_sha256
```

### Per-station section

Include:

``` text
station_id
tart_id
station_name
latitude_deg
longitude_deg
altitude_m
archive_query_status
archive_query_start
archive_query_end
raw_files
observation metadata
measurement_set metadata
```

### Synchronisation section

Include:

``` text
actual_common_start_utc
actual_common_end_utc
actual_common_duration_seconds
three_station_coverage_fraction
station_pair_time_offset_statistics
three_station_valid
```

## 24. Dry-Run Mode

Implement `--dry-run`.

It should:

1.  load Stage-1 results;
2.  select the window;
3.  resolve stations;
4.  calculate padded query intervals;
5.  query archive availability if possible;
6.  report what would be downloaded;
7.  perform no large downloads or conversions.

## 25. Logging

Use Python `logging`.

Log:

-   selected satellite/window;
-   station identifiers;
-   archive query interval;
-   archive results;
-   downloads;
-   hashes;
-   native-data validation;
-   `tart2ms` conversion;
-   synchronization;
-   output paths;
-   failures.

All public functions/classes require docstrings specifying purpose,
inputs, outputs, units, assumptions and failure behaviour.

## 26. Software Structure

Extend the current repository:

``` text
src/tart_gnss/
├── stage2/
│   ├── __init__.py
│   ├── selection.py
│   ├── archive.py
│   ├── download.py
│   ├── native.py
│   ├── synchronise.py
│   ├── measurement_set.py
│   ├── manifest.py
│   └── cli.py
```

Install an entry point such as:

``` text
retrieve-tart
```

### Responsibilities

`selection.py`: load/select/validate Stage-1 windows.

`archive.py`: TART archive abstraction and station mapping.

`download.py`: retrieval, integrity checking and caching.

`native.py`: inspect HDF5/JSON observations.

`synchronise.py`: coverage, common interval and timing diagnostics.

`measurement_set.py`: wrap `tart2ms` and validate MS output.

`manifest.py`: write the reproducible Stage-2 manifest.

`cli.py`: argument parsing and orchestration.

## 27. Tests

Tests must not depend entirely on live services. Use mocked archive
responses and small fixtures.

Required tests:

### Stage-1 selection

-   valid rank;
-   invalid rank;
-   valid NORAD;
-   ambiguous NORAD/window.

### Padding

Verify query start/end arithmetic.

### Missing station

Mock two available stations and one unavailable station; verify
`three_station_valid = false`.

### Temporal intersection

For:

``` text
A: 10:00 -- 10:20
B: 10:03 -- 10:18
C: 10:01 -- 10:25
```

verify:

``` text
common: 10:03 -- 10:18
```

### No intersection

Verify failure when the three station intervals have no common overlap.

### Integrity

Verify deterministic SHA-256 calculation.

### Manifest

Validate required fields and station count.

### Measurement Set wrapper

Test successful and failed `tart2ms` conversion using mocks/fixtures.

## 28. Error Handling

Use explicit exceptions such as:

``` text
Stage1SelectionError
ArchiveUnavailableError
ObservationNotFoundError
DownloadError
ObservationValidationError
NoCommonTimeError
MeasurementSetConversionError
```

Do not use bare `except:` and do not silently suppress
archive/conversion failures.

## 29. Terminal Summary

Successful runs should report:

``` text
TART GNSS — STAGE 2

Satellite:
  GSAT0219 (GALILEO 23)
  NORAD: 43566

Predicted common window:
  <start> -- <end>

Station                 Archive    Files    Data interval        MS
-------------------------------------------------------------------
Mauritius               OK         ...      ...                  OK
Madagascar              OK         ...      ...                  OK
South Africa            OK         ...      ...                  OK

Actual three-station overlap:
  start: ...
  end: ...
  duration: ...

Three-station valid: YES

Manifest:
  data/observations/.../manifest.json
```

If unsuccessful, state `Three-station valid: NO` and the specific
reason.

## 30. Stage-3 Interface

Stage 3 should require only the Stage-2 manifest:

``` bash
image-tart   --manifest data/observations/20261001/NORAD_43566/manifest.json
```

Stage 3 will independently image the three Measurement Sets and compare
detected source directions with Stage-1 predictions.

## 31. Do Not Implement Yet

Stage 2 MUST NOT implement:

-   WSClean imaging;
-   CLEAN optimisation;
-   image source detection;
-   satellite centroid fitting;
-   astrometric uncertainty estimation;
-   triangulation;
-   MCMC;
-   nested sampling;
-   Bayesian `(x,y,z)` inference;
-   TLE-vs-inferred residual analysis.

These belong to later stages.

## 32. Scientific Separation

Keep these concepts distinct:

**Geometric visibility** --- Stage 1 establishes that the satellite is
geometrically above the elevation threshold.

**Data availability** --- Stage 2 establishes that each TART actually
recorded usable visibility data.

**Radio detection** --- Stage 3 establishes whether the GNSS source is
detected in the interferometric image.

Do not infer one from another.

## 33. Initial Experimental Strategy

For the first end-to-end experiment:

1.  select one GPS or Galileo satellite;
2.  prefer high elevation at all three stations;
3.  choose a moderate-duration common window;
4.  check archive availability with `--dry-run`;
5.  initially retrieve only a short subsection;
6.  create three independent Measurement Sets;
7.  inspect timestamps, frequency and metadata;
8.  verify a genuine three-station temporal intersection;
9.  only then retrieve/process a longer interval.

## 34. Acceptance Criteria

Stage 2 is complete when:

``` bash
retrieve-tart   --visibility-results results/common_gnss_20261001.json   --rank 1   --padding 120   --output-dir data/observations   --download   --convert-ms
```

can:

1.  select exactly one Stage-1 satellite/window;
2.  resolve all three TART stations;
3.  query the appropriate archives;
4.  report whether observations exist at every station;
5.  retrieve native visibility files;
6.  preserve provenance and SHA-256 hashes;
7.  validate timestamps and observation metadata;
8.  determine the actual three-station temporal overlap;
9.  calculate Stage-1-window coverage;
10. independently convert each station to a Measurement Set;
11. validate each Measurement Set;
12. generate `manifest.json`;
13. explicitly identify partial/invalid three-station observations;
14. preserve sufficient provenance for reproducibility;
15. pass automated tests.

Scientific correctness, timestamp integrity, provenance and explicit
failure are more important than convenience or code brevity.
