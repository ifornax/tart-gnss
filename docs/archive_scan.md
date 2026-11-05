# Archive-aware candidate discovery (`tart-archive-scan`)

Stage 1 predicts **geometric** common visibility. Stage 2 checks the archive
for **one** selected window. A long Stage-1 search gives thousands of windows,
and many of them have no archived TART data at one or more stations. Checking
windows one by one with `tart-retrieve --dry-run` cannot show whether a
multi-station archive overlap exists *anywhere* in the interval.

`tart-archive-scan` works archive-first. It lists the archive once per
`(station, UTC day)`, turns the filenames into conservative per-station
*candidate* coverage, intersects the stations, and matches the shared
segments with **every** Stage-1 window. It is a utility, not a pipeline
stage. It never downloads a native file and never runs `tart-retrieve`.

Specification: [`tart_gnss_archive_scan_implementation.md`](../tart_gnss_archive_scan_implementation.md).

## 1. What a result means

These four evidence classes are kept separate. The scanner claims only the
second.

| Class | Evidence | Produced by |
|---|---|---|
| `GEOMETRICALLY_VISIBLE` | SGP4 prediction only | Stage 1 |
| `ARCHIVE_COVERAGE_CANDIDATE` | filenames and object metadata suggest coverage; no sample timestamp has been read | **this scanner** |
| `VALIDATED_THREE_STATION_DATA` | native HDF5 timestamps, frequency, identity and overlap validated | Stage 2 |
| `RADIO_DETECTED` | visibility-domain detection | Stage 3 (never inferred here) |

Every candidate row therefore carries:

- `evidence_class = ARCHIVE_COVERAGE_CANDIDATE`;
- `coverage_confidence = FILENAME_CANDIDATE`;
- `verified_overlap_duration_s` empty (null).

The `possible_overlap_duration_s` column is a planning estimate under an
explicit file-span assumption. It is **not** a measured integration time,
and it is not a proven upper bound.

## 2. Usage

Start with an inventory-only scan. It needs no Stage-1 input:

```bash
tart-archive-scan \
    --stations madagascar south_africa mauritius \
    --start 2026-10-09T00:00:00Z --end 2026-11-03T00:00:00Z \
    --min-common-duration 300 \
    --output results/archive_candidates.csv
```

Then match the archive with every Stage-1 window, reusing the cached listings:

```bash
tart-archive-scan \
    --visibility-results results/common_gnss_20261009.json \
    --min-common-duration 300 \
    --signal-relevance GPS_L1 GALILEO_E1 BEIDOU_B1C \
    --cache-policy reuse \
    --output results/archive_candidates.csv
```

| Option | Meaning | Default |
|---|---|---|
| `--stations ID [ID ...]` | station ids or aliases. With Stage-1 input, the set of canonical TART ids must equal the frozen Stage-1 stations, otherwise exit 2 | Stage-1 stations |
| `--start`, `--end` | aware UTC, half-open `[start, end)`. With Stage-1 input they may only **restrict** the Stage-1 interval | Stage-1 interval |
| `--visibility-results PATH` | Stage-1 JSON (or its CSV with the JSON alongside) | none: inventory-only |
| `--min-common-duration S` | minimum possible overlap of a reported candidate (inclusive) | 300 |
| `--signal-relevance V ...` | keep windows with these Stage-2 classifications | `GPS_L1` with Stage-1 input |
| `--constellations ...` | additional constellation filter | all |
| `--max-abs-tle-age-days D` | exclude candidates whose \|TLE age\| at the overlap midpoint exceeds D. The excluded count is reported; missing ages are flagged and kept | unset |
| `--station-file`, `--aliases` | catalogue and aliases. The catalogue supplies archive names only; coordinates come from Stage 1 | packaged |
| `--file-span-assumption-s S` | **assumed** maximum duration before a file's filename time | 120 |
| `--gap-tolerance-s S` | diagnostic grouping of segments in the manifest; gaps never count as covered | 0 |
| `--duration-cap-s S` | duration cap used in ranking | 600 |
| `--cache-dir`, `--cache-policy`, `--cache-max-age-hours` | inventory cache (`refresh`, `reuse`, `offline`) | `data/archive_inventory`, `refresh`, 24 |
| `--output PATH`, `--inventory-output PATH` | candidates CSV (companions are written alongside it), and the station/day CSV | `results/archive_candidates.csv` |
| `--min-stations N` | stations needed for a segment. Values below the number selected give **diagnostic** rows only | all selected |
| `--top N` | rows printed; every row is written | 20 |

**Exit status:**

- `0`: complete inventory, even if there are zero candidates.
- `3`: the scan finished, but the inventory is incomplete because a listing failed or a day is unknown.
- `2`: invalid arguments, station set, interval or Stage-1 input.
- `1`: an unexpected failure, or **every** listing failed.

**Exit classifications:**

- `CANDIDATES_FOUND`
- `COMMON_CANDIDATE_SEGMENTS_FOUND` (inventory-only mode)
- `NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS`
- `NO_IN_BAND_VISIBILITY_INTERSECTIONS`
- `ARCHIVE_LISTING_FAILURE`

When the inventory is incomplete, the classification gets an
`..._INVENTORY_INCOMPLETE` suffix, or becomes
`NO_CANDIDATE_OVERLAP_BUT_INCOMPLETE_INVENTORY`. The classification is
always based on segments that all selected stations cover; `--min-stations`
diagnostics never turn it into a "found".

## 3. Method

### 3.1 Listing (one call per station/day)

- **Reuse.** Station resolution, the archive backend, day prefixes and
  filename parsing all come from Stage 2. That means:
  - `resolve_run_stations` (look up `tart_id` first; identity and
    coordinates come from Stage 1);
  - `TartS3Archive`/`ArchiveBackend` (anonymous, read-only access, bucket
    `tart` at `s3.us-west-2.amazonaws.com`, location from `tart_tools`);
  - `day_prefixes` (`<archive>/vis/<Y>/<M>/<D>/`, no zero padding);
  - `parse_filename_utc`.

  Nothing is reimplemented.
- **Days.** A file named `T` can overlap `[start, end)` only if
  `start < T < end + span`. The scanner lists the days from `start − span`
  to `end + span`. This adds the day after `end` when `end + span` crosses
  midnight, and the day before `start` as an explicit edge allowance. For
  2026-10-09 → 2026-11-03 that is 27 days, so 81 prefixes for three
  stations.
- **Calls.** Each prefix gets exactly one successful `list_objects` call per
  scan. Transient errors are retried with bounded backoff (3 attempts, 2 s
  then 4 s). There are no requests per window and no downloads.
- **Canonicalisation.**
  - Only `vis_*.hdf|h5|hdf5` keys are kept; other keys are counted as
    ignored.
  - Duplicate keys are dropped.
  - Objects are sorted by `(filename_utc, key)`, so the output does not depend
    on listing order.
  - Size, ETag and S3 `last_modified` are kept, but `last_modified` (the
    upload time) is never used as observation time.
- **Status.**
  - `FOUND`: the listing succeeded and has at least one file.
  - `EMPTY`: the listing succeeded and has no files.
  - `ERROR`: the listing failed. An S3 error is never treated as an empty
    day.
  - `UNKNOWN`: the day was not listed (offline mode with no cache).
  - A station with no `archive_name` gives `ERROR` days.

### 3.2 Cache

- **Layout.** Snapshots are stored as
  `<cache>/<host>/<bucket>/<archive>/vis/<Y>/<M>/<D>/listing_<completed UTC>.json`.
  Each one records: `schema_version`, station, archive name, bucket, host,
  prefix, listing start and completion times, status, `objects` (key, size,
  ETag, last-modified, filename UTC), ignored keys, source, software version,
  and `checksum_of_canonical_listing` (SHA-256).
- **Writing.** Writes are atomic (a temporary file, then `os.replace`).
  `refresh` always writes a *new* snapshot and keeps the old ones as lineage.
  A failed listing is **never cached**.
- **Reading.** The schema and checksum are checked on read. A corrupt or
  incompatible newer snapshot is skipped with a warning and counted, and the
  next older valid snapshot is used.
- **`reuse`.** A snapshot is used only if both of these hold:
  - it is younger than `--cache-max-age-hours`;
  - it was made at least one hour after the day's UTC end (the day had
    "settled"; earlier listings can miss late uploads).

  Other days are relisted live.
- **`offline`.** Only the cache is used, the TTL is ignored, and a day
  without a valid snapshot is `UNKNOWN`. No S3 client is constructed.
- **Manifest.** The manifest records every snapshot's checksum and path, the
  number of live listings and list calls, and the cache hits, stale entries,
  corrupt files, failures and unknown days.

### 3.3 Candidate coverage, without filling gaps

- A file named `T` gives the candidate interval `[T − span, T]`. The default
  span is 120 s, the same as the Stage-2 end margin.
- Stage 2 found that the filename time is the last-sample time and that a
  file holds about 61 s of data. That observation is what motivates the
  span, but it is **not** a guarantee.
- A station's intervals are unioned only where they overlap or touch, so a
  gap is never covered. The result is clipped to `[start, end)`.
- Stations are intersected with Stage-2 `intersect_segments`, giving sorted,
  disjoint shared segments. Zero-length boundary contacts are discarded.
- **`--min-stations k < n`.** A sweep over all segment boundaries gives
  maximal intervals with a fixed station set, covered by at least `k`
  stations. These rows are labelled `TWO_STATION_DIAGNOSTIC` or
  `PARTIAL_NETWORK_DIAGNOSTIC` and are never valid retrieval inputs.
- **File counts.** A segment's count per station uses the Stage-2 selection
  rule: files with `start ≤ T ≤ end + span`.
- **Incomplete inventory.** A segment is `INVENTORY_INCOMPLETE` if any
  selected station's listing of a day that could hold files for it is not
  `FOUND` or `EMPTY`.

### 3.4 Matching with Stage 1

- **Loading.** The Stage-1 JSON is loaded with the Stage-2 loader and hashed
  once. Windows are grouped by NORAD ID and ordered by their `start_utc`
  string, exactly as `select_window` does, so `window_index` agrees with
  Stage 2 (a test checks this against `select_window`).
- **Filters.** Signal relevance comes from `stage2.selection.signal_relevance`:
  - GPS → `GPS_L1`;
  - Galileo → `GALILEO_E1`;
  - BeiDou-3 MEO/IGSO → `BEIDOU_B1C`;
  - BeiDou-2 → `OUT_OF_BAND`;
  - BeiDou-3 GEO and GLONASS → `UNKNOWN`, never treated as in band.

  Constellation filters are applied too. Every exclusion is counted in the
  manifest.
- **Matching.**
  - Each window is bisected into the sorted shared segments, so the cost is
    `O(W log S + K)`.
  - A long window that crosses several segments gives **one candidate per
    segment**; segments are never joined across a gap.
  - Observing-time totals come from the disjoint segments themselves, not
    from the candidate rows, so a segment seen by several satellites is
    counted once.
- **TLE age.** The age is `overlap midpoint − TLE epoch` in days (the
  Stage-1 `age_*_days` sign). The epoch is that satellite's own element-set
  epoch (`windows[].tle.epoch`), not the catalogue retrieval time. A missing
  epoch is flagged `TLE_AGE_MISSING`.

### 3.5 Ranking

The sort is deterministic. Each component is published per row in the JSON
and the rule is published in the manifest:

1. a complete inventory;
2. all selected stations present;
3. an in-band signal (`GPS_L1`, `GALILEO_E1`, `BEIDOU_B1C`; `UNKNOWN` is
   never preferred);
4. `min(possible duration, 600 s)`, descending. Because of the cap, a
   geostationary or IGSO window cannot win on length alone;
5. |TLE age| ascending;
6. Stage-1 geometry score, descending. This is the minimum Stage-1
   mid-window elevation over the stations;
7. Stage-1 rank, then NORAD ID, then overlap start, as tie-breakers.

### 3.6 Retrieval handoff

Each row has a ready-to-copy command for the **existing** Stage-2 CLI:

```bash
tart-retrieve --visibility-results <stage1.json> --rank <rank> --dry-run
```

`tart-retrieve --max-duration` centres its query on the *predicted window*
midpoint, not on the archive overlap, so it could miss a short overlap far
from the midpoint. The scanner therefore never emits `--max-duration`.
Each row reports the specific overlap interval to check against the dry-run
file listing instead.

## 4. Outputs

All outputs are written next to `--output`:

| File | Content |
|---|---|
| `archive_candidates.csv` / `.json` | ranked candidates (§3.4–3.6). The JSON also holds the run summary and the per-station counts as a map |
| `archive_candidate_stations.csv` | long form: one row per candidate × station, with the file count |
| `archive_common_segments.csv` | one row per shared segment (diagnostic rows included), with its inventory status and number of matched candidates |
| `archive_station_days.csv` | one row per station/day: status, files, bytes, earliest/latest filename, source, listing times, `covers_full_day`, checksum, cache path, error |
| `archive_station_summary.csv` | per station: days queried/found/empty/error/unknown, files, bytes, candidate intervals, possible covered seconds and fraction. An `ALL_STATIONS` row gives the intersection |
| `archive_scan_manifest.json` | inputs, git commit, S3 endpoint, frozen stations with registry metadata, assumptions, filters and exclusions, ranking rule, cache policy, listing statistics, snapshot hashes, Stage-1 SHA-256, summary and classification. A previous manifest is renamed, never overwritten |

The station-specific `files_<station>` columns are generated from the
selected stations. None are hard-coded.

## 5. Real-data result: 2026-10-09 → 2026-11-03

The examples are in [`docs/examples/archive_scan_real/`](examples/archive_scan_real/).
They were generated at commit `0a910ce` (clean tree) on 2026-11-03 at about
14:45 UTC. The Stage-1 input is `../sat_position/results/common_gnss_20261009.json`
(3,323 windows; its SHA-256 is in the manifests). The cache itself is
git-ignored.

**Listing.** The inventory-only run listed 81 prefixes live, with 81 list
calls and no failures. The matching runs used the cache for 78 prefixes and
relisted live the 3 prefixes for 2026-11-03, which had not yet settled.

| Station | Days queried | FOUND | EMPTY | ERROR | Files | Possible covered | Latest filename |
|---|---|---|---|---|---|---|---|
| mg-tana | 27 | 14 | 13 | 0 | 9,529 | 23.5 % | 2026-10-29T07:14:16Z |
| za-rhodes | 27 | 25 | 2 | 0 | 32,448 | 88.4 % | 2026-11-02T04:48:26Z |
| mu-udm | 27 | 0 | 27 | 0 | 0 | 0 % | — |

**Mauritius has no data.** All `mu-udm` listings succeeded and are empty.
A direct check of the archive confirms this: `mu-udm/vis/2026/` holds months
1–6 only, and the last file is `vis_2026-06-26_10_41_37.188062.hdf`. The
inventory is **complete**, so the result is
`NO_CANDIDATE_OVERLAP_WITH_COMPLETE_LISTINGS` (exit 0): no three-station
archive overlap exists anywhere in the interval. This is consistent with
the earlier spot checks of rank 8 (Rhodes only) and rank 590 (2026-11-02
13:46–19:08, after the last Rhodes file of the interval).

**Two-station diagnostic** (`--min-stations 2`; mg-tana + za-rhodes, **not**
a valid Stage-2 input):

- 69 shared segments, totalling 131.05 h.
- 1,086 in-band intersections of at least 300 s: 380 GPS L1, 358 Galileo E1
  and 348 BeiDou B1C. Another 476 were shorter than 300 s.
- About 94 h of the two-station time falls between 2026-10-24 and
  2026-10-29. Most of the rest is on 10-09 (10.8 h) and 10-17 (22.4 h).

Examples:

| Satellite | Stage-1 rank | Interval (UTC) | Duration | TLE age |
|---|---|---|---|---|
| GPS BIIF-9 (PRN 26), NORAD 40534 | 627 | 2026-10-28T21:52 – 10-29T03:10 | 5.3 h | −3.6 d |
| GSAT0213 (Galileo 17) | 1868 | 2026-10-25T12:28 – 15:04 | 2.6 h | −0.2 d |
| BeiDou-3 IGSO-3 (C40) | 2701 | 2026-10-28T23:08 – 10-29T02:06 | 3.0 h | −0.3 d |

NORAD 40534 is the satellite of the 2026-06-25 run.

**Next action for the real campaign.** A three-station campaign with this
triad needs one of these:

- `mu-udm` resumes archiving;
- an older interval when all three stations archived. The 2026-06-25
  Stage-2 run had all three, so an inventory-only scan over June 2026 is
  the next step;
- a different third station.

`za-rhodes` stops at 2026-11-02T04:48Z, and 11-03 was empty at listing
time. This could be an outage or upload lag; a later `--cache-policy
refresh` scan will tell which.

## 6. Validation

`tests/archive_scan/` has 54 tests, using a mocked `ArchiveBackend` that
records every listing and raises on download, and synthetic Stage-1 JSON.
No network access is needed.

The tests cover these areas:

- **Listing:**
  - prefix formatting and month/year transitions;
  - edge days;
  - exactly one listing per station/day over 25 days (81 calls);
  - parsing, malformed keys and duplicate keys;
  - the `FOUND`/`EMPTY`/`ERROR`/`UNKNOWN` distinctions;
  - bounded retry;
  - a missing `archive_name`.
- **Cache:**
  - refresh lineage;
  - reuse TTL;
  - relisting of a partial day;
  - offline mode with no S3 client;
  - checksum, schema and garbage rejection with fallback;
  - failures never cached;
  - separation of cache keys.
- **Intervals:**
  - gaps not filled;
  - touching merges;
  - half-open clipping across both midnights;
  - zero-length contacts;
  - the k-of-n sweep matching the intersection;
  - the file-count envelope;
  - gap groups.
- **Matching and ranking:**
  - no overlap despite geometric visibility;
  - an all-station overlap;
  - a two-station overlap that is not three-station valid;
  - separated segments kept separate;
  - the 300 s boundary (300 s kept, 299.9 s dropped);
  - sub-second partial overlap;
  - TLE age from the satellite epoch;
  - GPS L1 accepted, BeiDou-2 rejected, BeiDou-3 GEO `UNKNOWN`;
  - rank and window index equal to `select_window`;
  - long IGSO/GEO windows not dominating;
  - no double counting across neighbouring windows.
- **CLI:**
  - help and the public name;
  - alias resolution with frozen Stage-1 geometry;
  - station-set mismatch and invalid-input exits (2);
  - inventory-only mode;
  - a retrieval command that Stage 2's own parser accepts, without
    `--max-duration`;
  - `tart-retrieve` never run and no download;
  - an empty archive reported honestly;
  - a partial failure giving incompleteness and exit 3, including a
    candidate next to a failed day;
  - a total failure giving exit 1;
  - offline mode giving `UNKNOWN`;
  - output independent of listing order;
  - the previous manifest kept.

The tests also pass on Python 3.9.

## 7. Limitations

- Candidate coverage relies on the file-span assumption. Native files can
  contain internal gaps, which only Stage 2 can see.
- `verified_overlap_duration_s` is always null. Optional local native
  verification (`--verify-local-native`, spec §13) is **not implemented**.
- The handoff is a full-window dry run. Stage 2 has no
  `--query-start-utc`/`--query-end-utc` yet, so a short overlap far from
  the window midpoint has to be checked by hand in the dry-run listing.
- Long Stage-1 searches made with one current TLE catalogue propagate
  backwards for weeks. The TLE ages in the real run reach about −25 d. Use
  contemporaneous element sets for science.
- The one-hour "settled" margin used to judge whether a day's listing was
  complete is a heuristic for upload lag, not a guarantee.
- The edge-day allowance assumes a file's day directory is never *later*
  than the day after its filename time.

## 8. Issue, decisions and handoff

### Issue (ticket text)

Tracked as [issue #1](https://github.com/ifornax/tart-gnss/issues/1). Summary:

> **Archive-first candidate discovery (`tart-archive-scan`).**
>
> Stage-1 windows are geometric. Stage-2 dry runs check one window at a time,
> so the 3,323 windows from 2026-10-09 to 11-03 could not be ruled in or out.
>
> **Acceptance:**
> - one listing per station/day; no downloads; no `tart-retrieve` calls;
> - `refresh`/`reuse`/`offline` caching with checksums;
> - errors never reported as empty data, and an incomplete inventory gives
>   exit 3;
> - results labelled as candidates only;
> - Stage-1 rank, window index and NORAD ID preserved;
> - relevance and per-satellite TLE age reported;
> - frozen triad enforced;
> - deterministic ranking with a duration cap;
> - mocked-archive tests, and no regressions.
>
> **Failure modes:**
> - missing `archive_name`;
> - missing or mismatched Stage-1 input;
> - malformed or duplicate objects;
> - empty prefixes;
> - network errors;
> - a corrupt cache;
> - a mixed live/cache inventory;
> - partial failure;
> - invalid intervals, filters or policies.

### Decision log

1. **Filenames give candidates, not timestamps.** A filename is the
   last-sample time of files inspected so far, but it says nothing about
   gaps inside a file or about the true file length. Coverage built from
   filenames is therefore labelled `FILENAME_CANDIDATE`, and the verified
   duration stays null. Only Stage-2 native validation may promote it.
2. **The default span is 120 s, the Stage-2 margin, named as an assumption.**
   It is a search envelope. The documentation never calls it a hard bound.
3. **Prior-day allowance.** The scanner lists from `start − span`. This
   costs one extra prefix per station and protects against a file stored
   under an earlier day than expected.
4. **Failed listings are never cached.** This is the only way to guarantee
   that a cached day can never read back as `EMPTY` after an error.
5. **Partial days are stale in `reuse`.** A listing made before a day ended
   (plus one hour for upload lag) is relisted, even within the TTL.
6. **Stage-1 station resolution reuses `resolve_run_stations`.** The scanner
   behaves exactly as Stage 2 does, including legacy Stage-1 ids such as
   `mauritius`. Comparing canonical TART-id sets rejects a different triad.
7. **The classification is based on all-station segments only.**
   `--min-stations` diagnostics are reported separately and can never make
   the scan "found".
8. **The duration cap is 600 s.** Without it, GEO and IGSO windows lasting
   hours would always outrank MEO passes. With it, TLE age and geometry
   decide between candidates that are long enough.
9. **No `--max-duration` in the handoff.** That option centres on the
   predicted window, not on the overlap. Adding a Stage-2 interval option
   is left as a separate, backwards-compatible change.
10. **The scanner has no pandas dependency.** CSVs use the stdlib `csv`
    module, so the scanner runs in the core and `stage2` environments.

### Handoff

- **Branch and base.** Branch `archive-scan`, from `main` at `2b39795`
  (v0.8.0). Version is now 0.9.0.
- **Code.** The scanner is in `src/tart_gnss/archive_scan/`:
  - `inventory.py`, `cache.py`, `intervals.py`, `matching.py`,
    `ranking.py`, `output.py`, `errors.py`, `cli.py`;
  - entry point `tart-archive-scan` in `pyproject.toml` and
    `tests/test_cli_names.py`.
- **Stage 2 is unchanged.** Its semantics are untouched, and its public
  helpers are reused directly.
- **Tests.**
  - Run `pytest tests/archive_scan` (54 tests, about 3 s), then the full
    suite with `pytest`.
  - The full suite needs an environment with pyarrow (for the campaign
    Parquet tests) and minio. `test_real_tart2ms_conversion` needs a working
    `tart2ms` on `PATH`; this was a pre-existing environment failure.
- **Commands.** See §2. Use `--cache-policy reuse` after a first `refresh`.
- **Next actions.**
  1. Find an interval where all three stations archived. Run an
     inventory-only scan over 2026-06 first, or pick another third station.
  2. Optionally add Stage-2 `--query-start-utc`/`--query-end-utc`, checked
     to lie inside the Stage-1 window, so the scanner can emit commands for
     exact overlaps.
  3. Optionally add `--verify-local-native`, which would populate
     `verified_overlap_duration_s` as `TIMESTAMP_OVERLAP_VERIFIED` and fail
     closed.

## 9. Network-wide discovery (`--all-stations`)

Specification: [`tart_gnss_network_archive_discovery_implementation.md`](../tart_gnss_network_archive_discovery_implementation.md).
Code: `src/tart_gnss/archive_scan/network.py` (eligibility, sites,
baselines, sweep, triads) and `network_run.py` (run, tables, manifest,
report).

The explicit mode needs a triad chosen in advance. The network mode finds
**when and where the archive suggests that three or more geographically
distinct TART stations may have recorded simultaneously**, without any
Stage-1 input. It is still archive discovery between Stages 1 and 2, not a
new pipeline stage, and it reports `FILENAME_CANDIDATE` coverage only.

Four concepts stay separate:

1. a *registered* station (registry entry): not proof of operation or data;
2. a *filename-derived archive candidate*: what this mode reports;
3. *native-validated* simultaneous data: Stage 2;
4. a *detected* satellite: Stage 3.

### 9.1 Usage

```bash
tart-archive-scan --all-stations \
    --start 2026-07-01T00:00:00Z --end 2026-11-03T00:00:00Z \
    --min-stations 3 --min-common-duration 300 --min-baseline-km 500 \
    --cache-policy reuse --output results/network_overlap.csv
```

| Option | Meaning | Default |
|---|---|---|
| `--all-stations` | network mode; mutually exclusive with `--stations` | off |
| `--include-stations ID ...` / `--exclude-stations ID ...` | restrict / drop eligible stations (aliases accepted; unknown ids or an id in both lists → exit 2) | all eligible |
| `--min-stations N` | distinct physical **sites** per network interval, and the size of the reported station sets | 3 |
| `--min-common-duration S` | minimum **contiguous** candidate duration of a station-set segment | 300 |
| `--min-baseline-km KM` | every pairwise geodesic baseline of a set must be ≥ KM | no cut |
| `--site-tolerance-km KM` | stations within KM (geodesic, single linkage) are one site | 1 |
| `--max-triads N` | export/print only the top N sets; counts and scanning are unaffected | all |
| `--start`, `--end` | required; half-open `[start, end)` | — |
| `--file-span-assumption-s`, cache options, `--station-file`, `--aliases`, `--top` | as in the explicit mode | — |

`--visibility-results` (and the Stage-1 filters) are **rejected** in this
mode: a Stage-1 JSON is frozen to one station set and its geometry cannot be
applied to other triads. The network-only options are rejected without
`--all-stations`.

### 9.2 Method

1. **Eligible stations.** The catalogue is read entry by entry. Entries
   with invalid coordinates (including `metadata.skipped_entries` of a
   generated registry) are `INVALID_COORDINATES`; entries without an
   `archive_name` are `NO_ARCHIVE_NAME`; a repeated `archive_name` is
   `DUPLICATE_ARCHIVE_NAME` (listed once). `archive_verified` is recorded,
   never treated as evidence of data.
2. **Physical sites.** Union-find over all pairs with geodesic separation
   ≤ the tolerance. Group ids are `SITE_<lowest station id>`.
3. **Inventory.** The same `build_inventory` as the explicit mode: each
   eligible `(station, day)` prefix is listed (or read from the cache) once,
   with a progress line every 100 prefixes. The number of calls depends only
   on stations × days, never on the number of triads.
4. **Sweep.** One event sweep (`k_station_segments` with `k = 1`: half-open,
   ends before starts at equal instants, zero-length pieces dropped) gives
   maximal intervals with a constant active station set. Those whose active
   stations span ≥ N site groups are the **network active segments**
   (written even when a baseline cut leaves no triads). An active segment's
   `inventory_complete` considers **all** eligible stations, because a
   missing listing could hide a further station.
5. **Station sets.** Only sets of N stations at N distinct sites that are
   co-active in some active segment are considered. Each set's contiguous
   segments are the **intersection of its own members' segments**. A triad
   therefore stays one segment while a fourth station joins or leaves, and
   a co-located sibling never lends coverage (Stage 2 retrieves per
   station). The minimum duration applies to each contiguous segment;
   shorter pieces are counted (`n_short_segments_excluded`), never bridged.
6. **Baselines.** WGS84 Vincenty inverse geodesic (km), each pair computed
   once; the 3-D ECEF chord (with heights) is reported separately. The
   iteration raises on non-convergence (nearly antipodal points) rather
   than returning a wrong value.
7. **Ranking** (version 1, recorded in the manifest): complete inventory
   first; total possible shared duration; longest contiguous segment;
   minimum pairwise baseline (secondary only); station ids.

### 9.3 Outputs

For `--output results/network_overlap.csv`:

| File | Content |
|---|---|
| `network_overlap.csv` | ranked sets: durations, segment counts, first/last/longest segment, min/median/max baseline, chord, `coverage_class = FILENAME_CANDIDATE`, `verified_shared_duration_s` (null), proposed commands, `handoff_status = PROPOSED_NOT_RUN` |
| `network_overlap.json` | summary, the sets with their segments, and the baseline-rejected sets |
| `network_active_segments.csv` | sweep output: interval, active stations (`;` list and JSON), site groups, site count, completeness |
| `network_triad_segments.csv` | one row per set × contiguous segment: file counts per station, co-active additional stations, inventory status |
| `network_station_days.csv` | every station/day listing (status, files, source, checksum) |
| `network_station_summary.csv` | per station: coordinates, archive name, `archive_verified`, site group, day counts, files, first/last filename |
| `network_excluded_stations.csv` | excluded registry entries and reasons |
| `network_archive_scan_manifest.json` | arguments, software and git state, registry SHA-256 and metadata, archive endpoint, listing snapshots and checksums, cache statistics, interval, span assumption, site grouping, baseline method and pairs, policies, ranking rule, counts, output SHA-256s |

CSV and JSON outputs are replaced atomically; an earlier manifest is kept
under its run tag.

**Classifications:** `TRIADS_FOUND` (`_INVENTORY_INCOMPLETE`),
`THREE_SITE_OVERLAP_BUT_BASELINE_FILTER_EXCLUDES_ALL`,
`NO_THREE_SITE_CANDIDATE_OVERLAP_WITH_COMPLETE_INVENTORY`,
`OVERLAP_UNDETERMINED_DUE_TO_INCOMPLETE_INVENTORY`,
`NO_ELIGIBLE_ARCHIVE_STATIONS` (exit 0, empty tables) and
`ARCHIVE_LISTING_FAILURE`. "THREE" means `--min-stations` sites. More sites
requested than eligible sites exist → exit 2. Exit codes are those of the
explicit mode.

### 9.4 Handoff to the GNSS pipeline

The network scan never asserts GNSS visibility. For a chosen triad and
segment (`network_triad_segments.csv`):

1. run the proposed `tart-common-gnss --stations <triad> --start … --end …
   --min-elevation 20 --step 10 --tle-file <contemporaneous catalogue>`
   (the segment rounded outward to whole seconds);
2. run `tart-archive-scan --stations <triad> --visibility-results <new
   Stage-1 JSON>` to intersect in-band predictions with the archive;
3. `tart-retrieve --visibility-results <new JSON> --rank N --dry-run`;
4. download and validate (Stage 2) only viable candidates, then Stages 3–7.

Never reuse a Stage-1 JSON of another station set, and use element sets
contemporaneous with the segment.

### 9.5 Real-data result: 2026-07-01 → 2026-11-03

Run at commit `ddcaeed` (clean tree) with exactly:

```bash
tart-archive-scan --all-stations --start 2026-07-01T00:00:00Z --end 2026-11-03T00:00:00Z \
    --min-stations 3 --min-common-duration 300 --min-baseline-km 500 \
    --cache-policy reuse --cache-max-age-hours 240 \
    --output docs/examples/network_archive_scan_real/network_overlap.csv
```
 Outputs, terminal report and log:
[`docs/examples/network_archive_scan_real/`](examples/network_archive_scan_real).

- **Registry:** 18 entries; 14 eligible archive stations at 13 sites
  (`nz-nelson` + `nz-nelson-32` are one site). Excluded: `bw-biust`
  (`INVALID_COORDINATES`), `nz-rock2f`, `nz-test` and `tart3-test`
  (`NO_ARCHIVE_NAME`).
- **Inventory:** 1,778 station/day prefixes (14 × 127 days, including the
  two edge days), 1,700 listed live and 78 from the cache on the first run.
  There were no failures and no unknown days, so the inventory is
  **complete**. The rerun for the example used only cache hits.
- **Stations with no files in the interval:** `mu-udm`, `na-unam`,
  `nz-dunedin`, `se-kth` and `zm-cbu`. All of them are `archive_verified`
  in the registry.

| Station | Days with data | Files in interval | Possible coverage | Last filename UTC |
|---|---:|---:|---:|---|
| bd-iub | 125 | 169,245 | 96.4 % | 2026-11-02 04:47 |
| za-rhodes | 123 | 165,117 | 94.0 % | 2026-11-02 04:48 |
| ghana | 84 | 111,016 | 63.2 % | 2026-11-02 04:47 |
| nz-nelson | 58 | 77,277 | 43.1 % | 2026-11-02 04:47 |
| za-hammanskraal | 43 | 51,483 | 29.3 % | 2026-11-02 18:07 |
| mg-tana | 65 | 49,724 | 28.5 % | 2026-10-29 07:14 |
| nz-nelson-32 | 30 | 41,792 | 23.2 % | 2026-11-02 04:47 |
| nz-elec | 10 | 9,776 | 5.6 % | 2026-09-03 02:38 |
| it-ant | 1 | 198 | 0.1 % | 2026-09-28 13:42 |

(Exact values for every station are in `network_station_summary.csv`.)

- **Sweep:** 860 network active segments with ≥ 3 sites, totalling
  2,632.6 h of possible time with no double counting. Up to 6 sites are
  active at once.
- **Triads:** 44 sets of 3 distinct sites have a contiguous segment of at
  least 300 s. All 44 pass the 500 km all-pairs cut; the shortest pair is
  `za-hammanskraal`–`za-rhodes` at 896 km. All 44 have complete inventory.

| # | Triad | Total possible | Longest contiguous | Min baseline |
|---:|---|---:|---:|---:|
| 1 | bd-iub + ghana + za-rhodes | 1,749.8 h (107 segments) | 91.0 h, 09-17 13:55 → 09-21 08:57 | 5,168 km |
| 2 | bd-iub + nz-nelson + za-rhodes | 1,214.5 h | 125.6 h | 9,264 km |
| 3 | bd-iub + ghana + nz-nelson | 1,022.2 h | 158.3 h | 9,835 km |
| 4 | ghana + nz-nelson + za-rhodes | 954.3 h | 139.8 h | 5,168 km |
| 5 | bd-iub + mg-tana + za-rhodes | 760.4 h (232 segments) | 40.5 h, 07-28 23:51 → 07-30 16:18 | 2,633 km |

- **The previous triad.** `mg-tana` + `za-rhodes` + `mu-udm` has **no**
  candidate overlap, because `mu-udm` has no files after 2026-06-26.
  `mg-tana` + `za-rhodes` does pair with other third stations:
  - `bd-iub` (rank 5);
  - `ghana` (rank 13, 356.6 h);
  - `nz-nelson` (rank 15);
  - `za-hammanskraal` (rank 17).

These are filename candidates only. Ranking by archive duration says
nothing about GNSS geometry: `bd-iub` + `ghana` + `za-rhodes` spans Asia,
West Africa and South Africa, and satellites common to all three must be
established by Stage 1.

**Next actions** (proposed; none has been run):

1. Run `tart-common-gnss --stations bd-iub ghana za-rhodes --start
   2026-09-17T13:55:35Z --end 2026-09-21T08:57:47Z --min-elevation 20 --step
   10` with a TLE catalogue contemporaneous with September 2026. Do the same
   for `bd-iub mg-tana za-rhodes` on 2026-07-28T23:51:45Z → 07-30T16:18:51Z,
   which keeps the two stations already characterised in Stages 3–7.
2. Run the explicit `tart-archive-scan --stations … --visibility-results
   <new JSON>` on each result, then `tart-retrieve --dry-run`.
3. Promote a candidate only after Stage-2 native validation.

### 9.6 Validation

`tests/archive_scan/test_network.py` (29 tests, no network). The synthetic
five-station demo (docstring of the test module) has a co-located pair
(`cc`/`cc-2`, 110 m), an outage (`dd`), a station without an archive name
and one with invalid coordinates. Hand-computed truth that the tests check:

- network intervals with ≥ 3 sites: 10:30–11:50, 4 800 s in 7
  constant-membership pieces;
- `aa+bb+cc` is **one** 1 800 s segment although `dd` joins at 11:10;
  `aa+bb+cc-2` is 3 600 s; `cc` and `cc-2` never share a triad;
- the 240 s `dd` piece is excluded by the 300 s minimum and counted;
- a 500 km cut removes exactly the three sets containing both `aa` and `dd`
  (222 km); 1 000 km also removes the `bb–dd` (~890 km) sets;
- 15 listing calls for 5 stations × 3 days, none for excluded stations;
  cache reuse makes none; offline without a cache is a listing failure;
- a station/day error is `ERROR`, not an empty day, gives exit 3 and flags
  only the affected sets; no overlap with an incomplete inventory is
  `OVERLAP_UNDETERMINED…`;
- identical output bytes for different S3 listing orders;
- the Flinders Peak–Buninyong geodesic (54 972.271 m), the equatorial
  degree and the WGS84 quadrant; every packaged-registry pair converges;
- output SHA-256s, the proposed commands, and the CLI conflicts.

### 9.7 Limitations

- Everything in §7 applies (filename candidates, 120 s span assumption).
- `nz-dunedin` and `nz-elec` are 1.39 km apart, so at the default 1 km
  tolerance they are **separate** sites; use `--site-tolerance-km 2` (or a
  baseline cut) to merge them. `nz-nelson`/`nz-nelson-32` share coordinates
  and are always one site.
- A baseline cut is a coarse site-diversity filter. Crossing angles,
  elevations, clocks and calibration decide triangulation quality, and are
  evaluated after Stage 1 (Stage 4).
- S3 absence means no archived file, not that a telescope was not observing.
- Sets of size N = `--min-stations` only; there is no separate set-size
  option.

### 9.8 Issue, decisions and handoff

**Issue.** Tracked as [issue #2](https://github.com/ifornax/tart-gnss/issues/2). Summary:

> **Network-wide archive discovery (`tart-archive-scan --all-stations`).**
> The archive scan checks one hand-picked station set, and the chosen triad
> can have no overlap at all (`mu-udm` ends on 2026-06-26).
>
> **Scope:**
> - registry-derived eligible stations, with exclusions reported;
> - physical-site grouping;
> - one listing per station/day;
> - one sweep for intervals with ≥ N sites;
> - triads built from intersected per-station coverage;
> - an all-pairs geodesic baseline cut;
> - deterministic ranking and proposed Stage-1 commands.
>
> **Acceptance:** the §15 checklist of the specification.
>
> **Known limitations:**
> - filename candidates only;
> - a baseline cut is not a triangulation metric;
> - S3 absence does not prove that a telescope was not observing;
> - `nz-dunedin`/`nz-elec` (1.39 km apart) count as separate sites at
>   1 km.

**Decision log**

1. **Extend `tart-archive-scan`, do not add a command.** The inventory,
   cache, filename parser and outputs are shared; `--all-stations`
   dispatches to `network_run.run_network`, and the explicit mode's
   behaviour is unchanged (its 54 tests pass unchanged; its CSV/JSON writes
   are now atomic, with the same bytes).
2. **Sites, not ids.** `--min-stations` counts physical sites; co-located
   receivers would otherwise inflate the count. Single linkage at 1 km (the
   spec's example) is the default, independent of `--min-baseline-km`.
3. **One sweep, then per-set intersections.** The sweep finds where ≥ N
   sites are active; triads are enumerated only from co-active stations
   there, and each set's segments come from intersecting its own members'
   coverage. Stitching sweep pieces would have split a triad when a fourth
   station joins; unioning a site's stations would have credited coverage
   no single station has.
4. **Vincenty, no new dependency.** Neither `geographiclib` nor `pyproj` is
   a dependency. Vincenty is sub-millimetre on the reference line and fails
   loudly (rather than silently) near antipodes; no registry pair is close
   to antipodal.
5. **Reject `--visibility-results` in network mode.** The spec allows either
   rejection or identical-triad matching; rejection is the simplest
   guarantee that a Stage-1 JSON is never applied to another triad.
6. **Active-segment completeness uses all eligible stations**, set-segment
   completeness only the set's members. A missing listing elsewhere cannot
   remove a set's coverage, but it can hide a further active station.
7. **Exclusion of no-archive stations happens before listing.** Otherwise
   every run would carry `ERROR` days and exit 3.
8. **Atomic CSV/JSON writes** (both modes); the manifest keeps its earlier
   versions under their run tags and now records output SHA-256s in network
   mode.

**Handoff**

- Branch `network-archive-discovery` from `main` at `0d66eb0` (v0.9.0);
  version 0.10.0.
- New: `archive_scan/network.py`, `archive_scan/network_run.py`,
  `tests/archive_scan/test_network.py`. Changed: `cli.py` (options,
  dispatch), `inventory.py` (progress log), `output.py` (atomic writes).
- Tests: `pytest tests/archive_scan` (83 tests); full suite as in §8.
- Next actions are listed in §9.5.
- The real example was generated at `ddcaeed` (clean); the listing cache
  stays in `data/archive_inventory/` (git-ignored).
