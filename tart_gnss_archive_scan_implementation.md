# TART-GNSS Archive-Aware Candidate Discovery
## Implementation specification: `tart-archive-scan`

**Repository:** `ifornax/tart-gnss` (private)  
**Target branch:** branch from current `main`  
**Type:** utility/workflow enhancement, not Stage 8  
**Public CLI:** `tart-archive-scan`  
**Specification file:** `tart_gnss_archive_scan_implementation.md`

## 1. Purpose and problem statement

Stage 1 (`tart-common-gnss`) finds **geometric visibility**, not archived observations.
Stage 2 (`tart-retrieve --dry-run`) checks the archive for **one selected satellite window**.
A long search can return thousands of predicted windows, while many contain no archived
TART observations for the required station triad.

Concrete example from the current workflow:

- Station triad: `mg-tana`, `za-rhodes`, `mu-udm`.
- Stage-1 interval: `2026-10-09T00:00:00Z` to `2026-11-03T00:00:00Z`.
- Stage 1 found 3,323 common geometric visibility windows.
- Stage-2 dry run for rank 8 found 1,301 files at Rhodes, zero at the other two sites.
- Stage-2 dry run for rank 590 (GPS PRN 15; 2026-11-02 13:46–19:08 UTC)
  found zero files for all three sites.

These checks **do not** demonstrate that no three-station archive overlap exists
anywhere across the entire 25-day interval.

Implement a scanner that first inventories **actual archive keys** for the stations,
then intersects archived file-coverage **candidates** with the Stage-1 satellite
visibility windows, and ranks viable retrieval targets.

**Do not conflate**:

1. `GEOMETRICALLY_VISIBLE`: Stage-1 SGP4 prediction only.
2. `ARCHIVE_COVERAGE_CANDIDATE`: filename/object-metadata evidence of coverage;
   not yet verified from sample timestamps.
3. `VALIDATED_THREE_STATION_DATA`: native HDF5 timestamps, frequency, integrity
   and overlap verified by the existing Stage-2 validation process.
4. `RADIO_DETECTED`: Stage-3 visibility-domain detection; never inferred here.

The scanner's default output can only claim status 2, not 3 or 4.

## 2. Mandatory public command naming

All public commands use the existing `tart-*` namespace:

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
tart-campaign
tart-archive-scan         # new
```

Keep existing commands and semantics unchanged. Add the entry point to
`pyproject.toml`; update `tests/test_cli_names.py`, README, and docs.
Never rename previous commands or rewrite historical manifests.

## 3. Inspect and reuse the existing repository

Before implementation, inspect actual schemas and reuse these modules:

- `tart_gnss.stage2.archive.TartS3Archive`, `ArchiveBackend`,
  `ArchiveObject`, `ObservationRecord`, `parse_filename_utc`,
  `day_prefixes`, `find_observations`.
- `tart_gnss.stage2.selection.load_stage1_results`,
  `Stage1Window`, `signal_relevance`, `parse_utc`.
- `tart_gnss.stage2.synchronise.intersect_segments`,
  `total_seconds`, `coverage_fraction`.
- `tart_gnss.stations` station registry, aliases and canonical IDs.
- `tart_gnss.tle` formatting and SHA-256 utilities.
- Existing Stage-2 native-file validators for optional later validation.

**Important discovered implementation details:**

- The S3 archive uses the public `tart` bucket at the endpoint configured
  by `tart_tools`, with fallback `s3.us-west-2.amazonaws.com`.
- Day keys are `/<archive_name>/vis/<year>/<month>/<day>/`, with non-zero-padded
  year/month/day segments exactly as implemented in Stage 2 (do not guess).
- Filenames follow `vis_YYYY-MM-DD_HH_MM_SS.ffffff.hdf`.
- Filename time is the **last sample** in the native file in previously inspected
  examples. File duration is approximately one minute, **not guaranteed**.
- Stage-2 `find_observations` selects files by filename time, not S3
  `last_modified`, correctly avoiding upload-lag bias.
- Stage-2's default end margin is 120 seconds. Do not silently reinterpret
  it as a guaranteed file duration or exact coverage proof.

**Required architecture:** Extract or reuse a **shared listing/cache helper**
from Stage 2 if needed; do not duplicate bucket, parser, or station resolution
logic in an incompatible scanner-only implementation.

## 4. CLI contract

Provide:

```bash
tart-archive-scan \
    --stations madagascar south_africa mauritius \
    --start 2026-10-09T00:00:00Z \
    --end 2026-11-03T00:00:00Z \
    --visibility-results results/common_gnss_20261009.json \
    --min-common-duration 300 \
    --signal-relevance GPS_L1 \
    --output results/archive_candidates.csv
```

Arguments:

| Option | Meaning | Default |
|---|---|---|
| `--stations ID [ID ...]` | Station IDs or aliases; for Stage-1 matching, the set must match its frozen triad | From Stage-1 JSON when present |
| `--start`, `--end` | Explicit aware UTC interval; half-open `[start,end)` | Stage-1 interval when present |
| `--visibility-results PATH` | Existing Stage-1 detailed JSON or CSV with adjacent JSON | Optional; inventory-only mode without it |
| `--min-common-duration SECONDS` | Minimum intersection duration to report as a candidate | 300 |
| `--signal-relevance VALUE [VALUE ...]` | Filter satellite windows by Stage-2 classification | `GPS_L1` preferred when matching; no filter for inventory-only |
| `--constellations ...` | Optional additional constellation filter | all |
| `--max-abs-tle-age-days D` | Optional quality cutoff using per-satellite TLE epoch | unset (flag, never silently discard) |
| `--station-file`, `--aliases` | Reuse existing frozen/packaged registry and alias mapping | project defaults |
| `--file-span-assumption-s` | Maximum **assumed** prior duration of a native file, candidate discovery only | 120 |
| `--gap-tolerance-s` | Maximum gap between candidate segments allowed for a *diagnostic* grouping; gaps never count as covered | 0 |
| `--cache-dir PATH` | Persistent archive-inventory cache | `data/archive_inventory` |
| `--cache-policy` | `refresh`, `reuse`, `offline` | `refresh` |
| `--cache-max-age-hours H` | TTL for `reuse` mode | 24 |
| `--output PATH` | Ranked candidates CSV | `results/archive_candidates.csv` |
| `--inventory-output PATH` | Separate station/day summary CSV | alongside output |
| `--min-stations N` | Number of stations needed for an overlap | number selected, defaults to 3 for the triad |
| `--top N` | Limit printed rows; **not** scanned rows | 20 |
| `--verbose` | Detailed logs | off |

Prefer `argparse`, consistent with current code. Boolean options follow
existing repository CLI practices.

**Input policy:** If Stage-1 JSON is supplied, match station identities
to its frozen station set and station coordinates. Use the registry only for
archive names and API details. Reject mismatched station selections rather than
silently changing Stage-1 geometry.

`--start/--end` may restrict the Stage-1 search interval, never extend beyond
it unless explicitly documented as inventory-only.

## 5. Archive listing algorithm: only O(stations × days) remote listings

For a 25-day interval and three stations, list each required
`(station, UTC day)` S3 prefix **once per run** (plus explicit edge handling).

Do not query S3 separately for each of the 3,323 satellite windows.
Do not download native observations just to discover filenames.

Algorithm:

1. Resolve selected station IDs and archive names; record registry snapshot SHA.
2. Generate UTC day prefixes covering the interval, including any prior day
   needed if files ending after midnight could overlap the start boundary.
3. Use the existing anonymous MinIO/S3 archive backend to list all keys per
   `(station, day)`, once; canonicalise by full object key.
4. Keep only keys whose last component matches the established `vis_*.hdf`
   filename parser; record ignored/malformed objects.
5. Parse filename UTC; preserve key, size, ETag, S3 last-modified and day.
6. Sort and deduplicate by full key; do not double-count day-edge listings.
7. Write a station/day inventory with count, earliest/latest filename UTC,
   aggregate bytes, status (`FOUND`, `EMPTY`, `ERROR`, `UNKNOWN`), cache source
   and listing timestamp.
8. Preserve listing errors per station/day. **An S3 error is not empty data.**
9. Include station/day coverage in reporting even when zero objects were found.

Use a bounded streaming iterator if a prefix can contain many objects;
avoid loading the same prefix repeatedly.

## 6. Cache and reproducibility

Persist one versioned JSON (or JSONL) listing per station/day, with:

```text
schema_version
station_id
archive_name
bucket
host
prefix
listing_started_utc
listing_completed_utc
listing_status
objects: key, size_bytes, etag, last_modified_utc, filename_utc
source: live/cache
software_version
checksum_of_canonical_listing
```

Cache keys must distinguish bucket/host/station/archive prefix/day. Use atomic
writes (temporary file then replace) and verify schema/checksum when reading.

Policies:

- `refresh`: query live S3; persist a new snapshot, retaining previous history
  or enough lineage for reproducibility.
- `reuse`: reuse valid cache younger than TTL; refresh old/missing day caches.
- `offline`: use only available cache; explicitly label missing days as
  `UNKNOWN`, not `EMPTY`.
- Never transform a failed listing into an empty successful cache.

Record input Stage-1 JSON SHA-256 and every inventory snapshot hash in the
scanner manifest so results are reproducible from saved inventory.

## 7. Filename-based candidate coverage: conservative and explicit

**Critical limitation:** filename timestamps alone do not prove that any
intermediate seconds were sampled. Some files may contain gaps.

For each file ending at `t_file` create a possible candidate interval:

```text
[t_file - file_span_assumption_s, t_file]
```

For a sequence of files, union these intervals **without filling gaps**.
The assumed 120-second span is a search envelope, not a verified exposure
time. Keep two quantities distinct:

- `possible_overlap_duration_s`: interval overlap under the explicit
  file-span assumption; this is an **upper-bound-like planning estimate**,
  not measured integration time.
- `verified_overlap_duration_s`: only populated after actual native-file
  timestamp inspection, otherwise null.

Do not describe the assumed possible overlap as an exact, mathematical upper
bound unless the specified file-span assumption has been validated as a
hard bound for the archives under consideration.

Also compute conservative structural evidence:

- `files_per_station`;
- `candidate_intervals_per_station`;
- intersection of all station candidate segments;
- union of shared candidate intervals (do not double-count windows);
- archive errors and inventory completeness.

For `--min-stations 2`, label results explicitly **two-station diagnostics**,
not valid three-station retrieval inputs.

## 8. Multi-station intersection

Use the existing interval intersection utility after normalising sorted,
non-overlapping station candidate segments.

For every shared segment:

```text
start_utc
end_utc
possible_duration_s
station_ids
file_counts_per_station
inventory_status
coverage_confidence = FILENAME_CANDIDATE
```

Discard zero-duration boundary intersections. `--min-common-duration` filters
candidate **segment duration**, but does not claim actual validated coverage.

When inventory days are missing or errored, mark overlapping candidates as
`INVENTORY_INCOMPLETE`; do not claim definitive lack of common coverage over
the entire search range.

Compute a separate coverage summary for each station and the all-station
intersection.

## 9. Intersect with Stage-1 GNSS windows

When `--visibility-results` is given:

1. Reuse the Stage-1 JSON loader and parse all windows (not only one `--rank`).
2. Confirm the requested station set matches Stage-1 frozen stations.
3. Reuse `stage2.selection.signal_relevance(constellation, satellite_name)`.
4. Apply signal relevance/constellation filters and optional TLE-age quality
   restriction; record exclusion counts.
5. Intersect each Stage-1 visibility interval with shared archive candidate
   segments using sorted two-pointer/sweep-line matching, not quadratic
   `window × all file` scans.
6. For every intersection, calculate:
   - original Stage-1 rank, window index, satellite, NORAD, constellation;
   - signal relevance and its reason;
   - predicted-window UTC bounds and duration;
   - candidate-overlap start/end and possible duration;
   - station file counts within the overlap's lookup envelope;
   - geometry metrics directly from Stage-1 if available;
   - signed TLE age at overlap midpoint (from exact Stage-1 per-satellite
     TLE epoch), with missing/ambiguous ages flagged;
   - no scientific claims about detection or timing precision.
7. Preserve separate candidate segments of a long visibility window; do not
   collapse them into one fabricated continuous overlap.

For an archived station interval intersecting multiple satellite visibility
windows, do not double-count it in network observing-hours totals.

## 10. Candidate ranking and retrieval handoff

Ranking must be transparent and deterministic. Do not let very long
geostationary visibility windows dominate solely by duration, and do not
prefer `UNKNOWN` frequency relevance by default.

Recommended primary order:

1. complete inventory for all required stations;
2. validated in-band signal relevance (`GPS_L1`, `GALILEO_E1`, `BEIDOU_B1C`);
3. possible overlap duration descending **up to a configurable cap** (e.g. 600 s);
4. absolute TLE age ascending;
5. geometry score from documented existing Stage-1 fields if available;
6. Stage-1 rank, NORAD, UTC start as stable tie-breakers.

Publish each ranking component and the actual ranking rule in the manifest.

The output must contain a ready-to-copy retrieval command for the **existing**
`tart-retrieve` CLI, e.g.:

```bash
tart-retrieve \
  --visibility-results results/common_gnss_20261009.json \
  --rank 590 \
  --dry-run
```

**Important handoff caveat:** existing `tart-retrieve --max-duration` centres
the selection on the *geometric visibility window*, not the discovered archive
overlap. A short overlap far from that midpoint can be missed.

Therefore the scanner must:

- report the **specific archive-overlap interval** independently;
- never emit a `--max-duration` command claiming it targets that interval
  unless Stage 2 has acquired a genuine window-start/window-end selection API;
- in the first version, emit the original full-window `--dry-run` command and
  explain that the user may manually select an interval after Stage-2 support
  exists;
- optionally propose a **separate backwards-compatible Stage-2 enhancement**
  `--query-start-utc` and `--query-end-utc`, checked to lie inside the
  original Stage-1 visibility window; if implemented, add tests and manifest
  provenance, and allow scanner to emit precise short-interval commands.

No automatic download/convert operation is permitted from the scanner unless
explicitly authorised in a later workflow design.

## 11. Outputs

Recommended:

```text
results/
├── archive_candidates.csv
├── archive_candidates.json
├── archive_common_segments.csv
├── archive_station_days.csv
├── archive_station_summary.csv
└── archive_scan_manifest.json
data/archive_inventory/
└── ... station/day snapshots and hashes ...
```

`archive_candidates.csv` minimum columns:

```text
candidate_id
rank_stage1
norad_id
satellite_name
constellation
signal_relevance
window_index
predicted_start_utc
predicted_end_utc
overlap_start_utc
overlap_end_utc
possible_overlap_duration_s
verified_overlap_duration_s  # null unless validated
stations
n_stations
files_mg_tana
files_za_rhodes
files_mu_udm
tle_epoch_utc
tle_age_days
inventory_complete
coverage_status
ranking_score_or_order
retrieval_command
```

**Do not hard-code the initial triad** as schema columns or the only possible
stations: ideally station counts are a structured JSON map and a separate
long-form per-candidate-station table. If adding the human-facing triad columns
for convenience, dynamically generate them and preserve a general schema.

`archive_station_days.csv`: one row per station/day/status/count/bytes.
`archive_station_summary.csv`: station, days queried, days FOUND/EMPTY/ERROR,
files, possible covered seconds, coverage %, earliest/latest filename.
`archive_common_segments.csv`: one row per candidate shared archive interval.

Manifest must contain all inputs, code/git version, S3 endpoint/bucket,
stations with frozen coordinates and registry metadata, assumptions, filters,
cache policy, listing status, hashes, candidate counts and exit classification.

## 12. Human-readable terminal report

Example shape (numbers are illustrative, never hardcode):

```text
TART ARCHIVE SCAN
Interval: ...
Stations: mg-tana, za-rhodes, mu-udm

Station    Days queried  Days with data  Files  Listing errors
mg-tana    ...           ...             ...    ...
za-rhodes  ...           ...             ...    ...
mu-udm     ...           ...             ...    ...

Possible three-station candidate segments: ...
Candidate segments >= 300 s: ...
In-band satellite intersections: ...
Inventory completeness: COMPLETE / INCOMPLETE
Coverage class: FILENAME_CANDIDATE (NOT NATIVE VALIDATED)

Top candidates:
  rank   NORAD   relevance   overlap start/end   possible duration   TLE age
  ...
```

When no overlap is found distinguish:

- `NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS`
- `NO_CANDIDATE_OVERLAP_BUT_INCOMPLETE_INVENTORY`
- `ARCHIVE_LISTING_FAILURE`
- `NO_IN_BAND_VISIBILITY_INTERSECTIONS`

Never print an unqualified “no common data” when inventory is incomplete.

## 13. Optional native-data verification (later / opt-in)

Default scanner does **not** download HDF5.

If `--verify-local-native PATH` is implemented, it may read already downloaded
files to verify:

- HDF5 format/shape and observation metadata;
- actual timestamps and data gaps (reuse Stage-2 native utilities);
- frequency includes 1575.42 MHz;
- valid station identity and coordinates;
- true three-station timestamp-segment overlap and cadence.

Only then set `verified_overlap_duration_s`. Do not call
`VALIDATED_THREE_STATION_DATA` unless **all** existing Stage-2 validation
requirements are satisfied. Otherwise use `TIMESTAMP_OVERLAP_VERIFIED`,
which is a narrower, accurate label.

Do not estimate actual coverage from filename spacing and call it validated.

## 14. Numerical/time correctness

- All internal timestamps timezone-aware UTC, including `Z`/offset input.
- Half-open query ranges, documented boundary inclusions.
- Include preceding day boundary when necessary.
- Sort and deduplicate all intervals before intersection.
- Avoid rounding 0.1 s Stage-1 window boundaries to whole seconds.
- TLE age must correspond to the *individual satellite element set*, not
  aggregate catalogue retrieval time.
- Zero-duration intersections are excluded.
- A file ending just after the query start can still contain earlier samples:
  include it in candidate discovery.
- Do not count S3 upload time as observation time.
- All coverage fractions have denominators defined and are in [0,1].
- Do not double-count overlapping station/file/window intervals.

## 15. Performance requirements

For 25 days, three stations and 3,323 predicted windows:

- one live listing per required station/day prefix per scan (plus documented
  edge cases);
- no requests per satellite window;
- output sorted deterministically;
- cache reuse avoids redundant S3 scans;
- interval matching approximately O((W + S) log(W + S)) or better after
  sorting, where W is visibility windows and S shared candidate segments;
- memory usage proportional to archive objects for the selected interval,
  not observation-sample arrays.

Log the number of prefixes scanned and cache hits.

## 16. Security, reliability and errors

Use anonymous read-only public S3 access and bounded retry/backoff for
transient errors. Do not require credentials. Never mutate the archive.

Handle:

- station missing `archive_name`;
- missing Stage-1 JSON;
- different station sets;
- malformed filename;
- malformed/duplicate objects;
- empty prefixes;
- network timeout/permission errors;
- corrupt/incompatible cache;
- mixed live/cache inventory;
- partial failure in one or more station-days;
- invalid date interval, filter or cache policy.

Use `logging.info` for major steps; debug for prefix-level details. Avoid
flooding logs with individual objects.

Suggested exit codes:

```text
0  complete successful scan (even if zero genuine archive candidates)
3  scan completed but inventory incomplete or candidates scientifically
   indeterminate; outputs written with warning
2  invalid arguments/input/schema
1  unexpected technical failure
```

Do not return exit `0` if there were unreported archive errors.

## 17. Tests (mandatory)

Create `tests/archive_scan/` with mocked `ArchiveBackend` and deterministic
synthetic Stage-1 JSON.

Minimum tests:

1. Station aliases and frozen-geometry consistency.
2. Exactly one listing call per station/day prefix (including 25-day test).
3. Day prefix formatting for one/two-digit months and year/month transitions.
4. Cross-midnight file ending after start included.
5. Filename parsing, malformed objects and full-key deduplication.
6. `EMPTY` vs `ERROR` vs `UNKNOWN` distinctions.
7. Cache `refresh` / `reuse` / `offline`, TTL, checksum validation.
8. No overlap despite geometric visibility.
9. All-station common candidate overlap.
10. Two-station-only overlap is not three-station-valid.
11. Multiple separated segments remain separated; gaps not counted.
12. Minimum common duration boundary at 300 s.
13. Search-window intersection with partial overlap and subsecond precision.
14. TLE age computed from satellite epoch, not catalogue time.
15. GPS L1 accepted, BeiDou-2 out-of-band rejected by filter,
    BeiDou-3 GEO unknown not silently treated as B1C.
16. Rank and window-index preserved from Stage 1.
17. Long GEO windows not allowed to dominate ranking purely by length.
18. No duplicate counting across neighbouring visibility windows.
19. Candidate retrieval command targets a valid existing Stage-2 CLI, and
    does not misuse centred `--max-duration`.
20. Scanner never calls `tart-retrieve` or performs native downloads.
21. Empty complete archive yields honest “no candidate” reporting.
22. Partial failed day yields incompleteness flag and nonzero exit.
23. CLI options/help and public `tart-*` naming.
24. Output deterministic independent of S3 listing order.
25. Any optional local-native verification tests include true gaps and
    fail closed when HDF5 validation is insufficient.

Run existing Stages 1–7 and campaign tests to ensure no regression.
Unit tests must not require network access.

## 18. Real-data acceptance protocol

First run an **inventory-only** scan over the recent interval:

```bash
tart-archive-scan \
  --stations madagascar south_africa mauritius \
  --start 2026-10-09T00:00:00Z \
  --end 2026-11-03T00:00:00Z \
  --min-common-duration 300 \
  --output results/archive_candidates.csv
```

Examine per-station archive-day availability **before** concluding there are
no common data.

Then intersect with the existing Stage-1 results:

```bash
tart-archive-scan \
  --stations madagascar south_africa mauritius \
  --start 2026-10-09T00:00:00Z \
  --end 2026-11-03T00:00:00Z \
  --visibility-results results/common_gnss_20261009.json \
  --min-common-duration 300 \
  --signal-relevance GPS_L1 GALILEO_E1 BEIDOU_B1C \
  --cache-policy reuse \
  --output results/archive_candidates.csv
```

The scanner should return honest results even when no common archive window is
found. Save the report and use actual source status rather than guessing.

For a top candidate, manually run the emitted `tart-retrieve --dry-run`,
then `tart-retrieve --download --convert-ms` only after confirming candidate
coverage and resource requirements.

Note: Stage-1 25-day geometric searches with one current TLE catalogue may
include stale backward propagation. Preserve the age warnings, and use
contemporaneous element sets for scientific validation.

## 19. Documentation and repository handoff

Create or update:

```text
src/tart_gnss/archive_scan/
    __init__.py
    cli.py
    inventory.py
    cache.py
    intervals.py
    matching.py
    ranking.py
    output.py
    errors.py
tests/archive_scan/
docs/archive_scan.md
README.md
pyproject.toml
tests/test_cli_names.py
```

Adjust module structure to existing conventions without gratuitous duplication.

**Repository engineering requirements:**

- Create an issue/ticket describing the missing archive-first workflow,
  acceptance criteria, and failure modes.
- Maintain implementation notes and a decision/history log (including why
  filenames produce *candidates*, not verified timestamps).
- Update agent handoff documentation with current status, commands, tests,
  known limitations and next actions.
- Record exact base commit and branch used.
- Every public function must have clear docstrings specifying purpose,
  parameters, outputs and logic; use logging appropriately.

Do not change any existing Stage-1/2 scientific semantics without explicit
backwards-compatible specification and tests.

## 20. Final acceptance checklist

- [ ] `tart-archive-scan` installed and in CLI naming test.
- [ ] Inventory-only scan works without Stage-1 JSON.
- [ ] Stage-1 all-window matching works.
- [ ] Stations resolve to canonical TART identifiers and frozen geometry.
- [ ] S3 archive listing is reused rather than reimplemented incompatibly.
- [ ] No file is downloaded during default scan.
- [ ] At most one listing per unique station/day in a scan.
- [ ] Station/day cache supports refresh/reuse/offline and provenance.
- [ ] Archive listing errors never become `0 files` with a false success.
- [ ] Per-station availability and possible shared segments are reported.
- [ ] Signal relevance and per-satellite TLE age are included.
- [ ] Ranked candidates preserve Stage-1 rank/NORAD/window index.
- [ ] Retrieval commands are safe and compatible with Stage 2.
- [ ] Missing data and no candidate overlap are correctly distinguished.
- [ ] Candidate coverage is never misreported as validated observation data.
- [ ] Tests include boundary, gaps, stale-cache and failure cases.
- [ ] All existing tests still pass.
- [ ] README, docs, issue log, decision log and handoff are updated.

### Final agent report

Report:

```text
Git commit / branch
Files added and modified
CLI help and example commands
Archive prefixes listed / cache hits / failures
Station/day availability table
Number of possible shared segments
Number of intersected GNSS candidate windows
Top GPS L1 / Galileo E1 / BeiDou B1C candidates
Inventory completeness
Synthetic tests passed/failed/skipped
Full test suite passed/failed/skipped
Known limitations
Next action for the real TART campaign
```

**Scientific acceptance:** The scanner is complete when it can determine
where there are plausible multi-station archived observing intervals, and
match them to Stage-1 satellite predictions, **without falsely claiming
verified observation overlap or radio detection**.
