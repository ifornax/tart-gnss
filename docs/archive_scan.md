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
