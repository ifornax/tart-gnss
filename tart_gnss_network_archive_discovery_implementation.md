# TART-GNSS — Network-Wide Archive Discovery Implementation

**Deliverable:** Extend the existing `tart-archive-scan` command with `--all-stations`, configurable station-count and geographic-baseline constraints, and network-wide archive overlap discovery.

**Repository:** `ifornax/tart-gnss` (private)  
**Target branch:** create a feature branch from current `main`  
**Classification:** Archive discovery capability between Stages 1 and 2; **not a new pipeline stage**  
**Public CLI:** `tart-archive-scan` (unchanged)  
**Proposed specification filename:** `tart_gnss_network_archive_discovery_implementation.md`

## 1. Objective

Find time intervals during which **at least three geographically distinct TART stations may have archived, simultaneous visibility observations**. Do so efficiently across the complete station registry, without first selecting a satellite or enumerating all station triads and independently querying the archive for each triad.

Output candidate multi-station time segments and ranked triads that can subsequently be tested using `tart-common-gnss` and validated through `tart-retrieve`.

Keep four concepts separate:

1. **Installed / registered station**: an entry in the TART registry. This does **not** prove operation or archive availability.
2. **Filename-derived archive candidate**: S3 keys suggest an observation interval. This does **not** establish continuous recording or actual timestamp overlap.
3. **Native-validated simultaneous observation**: measured HDF5 sample timestamps and metadata pass Stage-2 checks at three stations.
4. **Detected satellite**: Stage-3 fitting confirms a satellite signal. Neither archive scan nor Stage 1 establishes this.

The initial implementation is **metadata-only**; do not download HDF5 during discovery, and never label filename-derived overlaps as verified observations.

## 2. Existing repository integration — inspect first

Before modifying code, read the actual implementations and schemas on the current branch, particularly:

- `src/tart_gnss/archive_scan/` — CLI, inventory, cache, intervals, matching, ranking, output and error handling (confirm filenames before editing).
- `src/tart_gnss/stage2/archive.py` — `TartS3Archive`, `ArchiveBackend`, `ArchiveObject`, `ObservationRecord`, `parse_filename_utc`, `day_prefixes`.
- `src/tart_gnss/stage2/synchronise.py` — interval intersection/coverage utilities.
- `src/tart_gnss/stage2/selection.py` — Stage-1 JSON, signal relevance, timestamp parsing.
- `src/tart_gnss/stations.py` and `src/tart_gnss/config/stations_registry.yaml` — station records, archive identifiers, coordinates.
- `src/tart_gnss/config/station_aliases.yaml` — friendly CLI aliases.
- `pyproject.toml`, `README.md`, `docs/archive_scan.md`, `tests/archive_scan/`, `tests/test_cli_names.py`.

Current frozen registry (retrieved 2026-10-09) contains 17 stations, but some have `archive_name: null`; registry presence, `archive_verified: true`, and API URL presence do not guarantee HDF5 data in the requested period. The registry may evolve: **derive the station set dynamically**, do not hardcode these counts.

Use the same anonymous read-only S3 bucket, archive-prefix convention, filename parser, snapshot/cache provenance and logging as the existing scanner. Maintain all existing `tart-*` commands and ensure the old explicit-station mode is backward compatible.

## 3. Primary proposed CLI

```bash
tart-archive-scan \
  --all-stations \
  --start 2026-07-01T00:00:00Z \
  --end 2026-11-03T00:00:00Z \
  --min-stations 3 \
  --min-common-duration 300 \
  --min-baseline-km 500 \
  --cache-policy reuse \
  --output results/network_overlap.csv
```

**These options are an implementation target, not a claim that the current CLI supports them.**

### CLI behavior

| Argument | Behavior | Default |
|---|---|---|
| `--all-stations` | Enumerate registry stations with usable `archive_name`; no explicit station list required | false |
| `--stations ID [ID ...]` | Preserve existing explicit-station workflow | unchanged |
| `--min-stations N` | Minimum **distinct geographic station sites** in a candidate shared interval | 3 in all-stations mode; retain existing defaults elsewhere |
| `--min-common-duration SECONDS` | Minimum contiguous *candidate* shared-segment duration | 300 |
| `--min-baseline-km KM` | Minimum **pairwise geodesic baseline** among every member of a reported triad, unless a separately named policy is selected | unset (no geographic cutoff) |
| `--max-triads N` | Optional limit on **displayed/exported ranked triads**, after network discovery; do not limit archive scanning | unset |
| `--include-stations IDs ...` | Restrict eligible registry stations to provided subset in all-stations mode | all eligible |
| `--exclude-stations IDs ...` | Exclude certain sites | none |
| `--station-file`, `--aliases` | Existing registry/alias behavior | existing |
| `--start`, `--end` | Time interval as aware UTC timestamps, half-open `[start,end)` | required for all-stations mode |
| `--visibility-results PATH` | Optional Stage-1 matching, only when frozen station set exactly matches an emitted candidate triad | none |
| `--file-span-assumption-s` | Existing assumed possible interval before each filename timestamp | existing default, currently 120 s |
| `--gap-tolerance-s` | Existing diagnostic grouping semantics; never convert missing spans into covered seconds | existing |
| `--cache-policy`, `--cache-dir`, `--cache-max-age-hours` | Reuse existing inventory cache | existing |
| `--output PATH` | Main CSV; write companion tables and JSON/manifest | `results/network_overlap.csv` |
| `--top N` | Number printed to terminal | existing |

`--all-stations` and `--stations` must be mutually exclusive unless a deliberate, tested semantics is chosen. `--include-stations` and `--exclude-stations` work only with `--all-stations` and should error on contradictions or unknown IDs.

**Important:** a Stage-1 file for `mg-tana+za-rhodes+mu-udm` cannot be applied to arbitrary triads; the frozen Stage-1 coordinates and geometry belong to that station set. In first release, `--visibility-results` in all-stations mode should either be rejected with an explanation, or explicitly match **only** the identical triad and output unmatched triads as archive-only. Never recompute Stage-1 geometry silently or use a mismatched prediction.

## 4. Eligible station discovery

1. Load the versioned station registry and aliases, preserving its source URI, retrieval time and SHA-256.
2. Include stations only if a usable `archive_name` is configured and coordinates are valid. `archive_verified` is a configuration assertion, not evidence of current archive availability.
3. Report excluded registry entries and reasons (`NO_ARCHIVE_NAME`, `INVALID_COORDINATES`, etc.).
4. Group exact duplicate/near-identical physical locations carefully: multiple receivers/testbench IDs at the same coordinates must not count as geographically independent sites. Default site grouping should use a documented small separation tolerance (e.g. 1 km), independent of `--min-baseline-km`; retain each station's identity and report grouping details. Configurable threshold optional.
5. Do not confuse the South African `za-rhodes` and `za-hammanskraal` stations; both are distinct.
6. Identify geodetically separate observatories before triad generation.

### Station summary table

Include station ID, name, latitude, longitude, altitude, archive name, `archive_verified`, date coverage, number of objects, earliest/latest filename UTC in scan interval, inventory completeness, number of days queried/with files/errors/unknown, and site-group ID.

## 5. Inventory once per station/day, not once per triad

Reuse the existing `TartS3Archive` `list_objects()` abstraction and `archive_scan` inventory/cache utilities. For each eligible station:

- Query the relevant UTC day prefixes (including edge prefixes required by the file-span assumption).
- List each unique `(bucket, prefix)` at most once during a run; use the checksummed cache according to the existing policy.
- Parse the filename-derived observation end timestamps; do not use S3 `last_modified` as observing time.
- Distinguish `FOUND`, `EMPTY`, `ERROR`, `UNKNOWN` and cached/stale states.
- Deduplicate by full S3 object key and preserve metadata (`size`, `etag`, `last_modified_utc`, timestamp).
- An S3 listing failure must **not** appear as a legitimate zero-file day.
- Do not download HDF5 in default network discovery.

A 125-day period with 14 eligible station prefixes is approximately 1,750 prefix/day listings before cache reuse; implement bounded I/O and useful progress logs rather than needlessly listing for every station triad.

The exact number of prefixes must come from the registry and boundary handling, not a hardcoded estimate.

## 6. Convert archive filenames to *possible* coverage

Maintain the existing semantics of `tart-archive-scan`: each filename timestamp `t` denotes a candidate file ending at `t` and provides an **assumed possible coverage** interval

\[
I_i(t)=[t-\Delta_{\mathrm{assumed}}, t],
\]

where `Δ_assumed = --file-span-assumption-s` (currently 120 s by default). Merge overlapping intervals for the same station, preserving genuine candidate gaps. Clip interval endpoints to the query range.

These are **candidate envelopes**, not validated HDF5 timestamp segments; avoid claims of exact exposure duration. If the assumption is wrong, a real segment can be missed or a false overlap generated. Record the assumption in every output and warning.

Report both `possible_duration_s` and `verified_duration_s` (null unless an optional true native-validation process is implemented).

## 7. Global multi-station overlap algorithm — required

Do **not** enumerate every possible triad first and separately intersect all of their intervals. Perform one **sweep-line / event-based network overlap** over the union of per-station candidate segments:

1. Emit a `+station` event at each start, and a `-station` event at each end for half-open intervals.
2. Sort all distinct UTC event instants, using deterministic tie-breaking for coincident events and half-open interval semantics; avoid zero-length spans.
3. Maintain the active **site groups** (not only station IDs).
4. Between consecutive instants, if `active_site_groups >= min_stations`, emit a network segment with start/end, duration, active stations and distinct-site counts.
5. Merge adjacent segments only if their active station/site sets and status are identical. Never bridge a real gap.
6. The minimum-duration filter applies to **contiguous intervals for a particular eligible triad**, not merely to brief global network segments that can change membership while a triad remains continuously present.
7. Keep incomplete-inventory periods flagged (or excluded from definitive zero-coverage findings). Do not report absence of coverage as verified if any required station/day listing is unknown.
8. Persist `network_active_segments` even where a minimum baseline filter leaves no eligible triads.

Suggested complexity: O(E log E) for E interval-boundary events, followed by triad enumeration **only within periods that actually have >=3 active distinct sites**.

## 8. Identify and rank geographically distinct triads

For each qualifying active station set, identify possible triads (or N-station sets if requested later). Use cached geodesic pairwise distances computed once per pair.

**Baseline definition:** geodesic/ellipsoidal surface distance between stations (WGS84) in kilometres; report separately the 3-D ECEF chord distance if helpful. Do not mix geodesic and chord units or methods.

For each triad:

- calculate all three pairwise baseline lengths;
- `min_pairwise_baseline_km`, `median_pairwise_baseline_km`, `max_pairwise_baseline_km`;
- require **all three pairwise** baselines >= `--min-baseline-km` when enabled;
- use exactly three distinct site groups; co-located stations do not qualify as independent triangulation anchors;
- merge consecutive time segments in which all three chosen stations are active, even if additional stations appear or disappear;
- compute individual and total possible shared duration without double-counting overlapping segments;
- attach station-archive file counts and coverage status;
- preserve the underlying time segments and their station IDs.

Do not imply that a 500 km minimum baseline guarantees good triangulation: satellite elevation, 3-D line-of-sight crossing angles, clock quality and calibration matter. Geometry is evaluated later using Stage 1/4, not at archive-only discovery time.

### Ranking

Create a transparent and deterministic ordering, for example:

1. `inventory_complete` true first;
2. total possible shared duration (descending);
3. longest continuous possible shared interval (descending);
4. minimum pairwise baseline (descending; **secondary**, not sole geometry proxy);
5. station IDs (lexicographic stable tie-break).

Allow explicit ranking components to be revised as scientific requirements emerge, but record the chosen rule/version in the manifest. Avoid a black-box score.

With 10+ simultaneously active stations, the number of potential triads may become large; group identical triads across time and write long-form interval records separately. `--max-triads` should limit export/printing, not change completeness statistics or hide required inventory work.

## 9. Recommended outputs

Given `--output results/network_overlap.csv`, create:

```text
results/
├── network_overlap.csv                 # ranked distinct triads
├── network_overlap.json                # machine-readable triads + intervals
├── network_active_segments.csv         # network-level event sweep output
├── network_triad_segments.csv          # one row per triad x continuous candidate interval
├── network_station_days.csv            # inventory per station/day
├── network_station_summary.csv         # availability per station
├── network_excluded_stations.csv        # excluded registry entries with reasons
└── network_archive_scan_manifest.json  # full provenance and quality controls
```

Minimum `network_overlap.csv` columns:

```text
rank
triad_id
station_ids
site_group_ids
n_stations
possible_shared_duration_s
longest_contiguous_duration_s
n_contiguous_segments
first_shared_start_utc
last_shared_end_utc
min_pairwise_baseline_km
median_pairwise_baseline_km
max_pairwise_baseline_km
inventory_complete
coverage_class                 # FILENAME_CANDIDATE
verified_shared_duration_s     # null until native validation
```

In `network_triad_segments.csv` store interval start/end/duration, triad, candidate confidence, source station file counts, inventory quality, and any co-active additional stations.

In `network_active_segments.csv` store active station IDs and distinct-site count for every nonzero candidate interval; use a JSON string or well-documented list encoding in CSV.

In the manifest preserve run ID, exact input arguments, software version, git commit/dirty state, registry SHA/snapshot, S3 endpoint/bucket, archive listing hashes and statuses, cache hits/age, UTC interval, file-span assumption, min duration, site grouping tolerance, pairwise baseline method, quality policies, ranking rule, counts and output hashes.

Do not silently overwrite incompatible previous scans: atomically replace explicitly requested outputs or create a versioned run subdirectory, consistently with existing pipeline conventions.

## 10. Terminal report

Display:

```text
TART NETWORK ARCHIVE DISCOVERY
Interval: ...
Registry stations: ...
Eligible archive stations: ...
Distinct geographic sites: ...
S3 station/day prefixes listed: ...  cache hits: ...  failures: ...

Station        Days with data  Files  First filename UTC  Last filename UTC
...

Network active intervals with >= 3 sites: ...
Ranked eligible triads: ...
Top triad          Continuous possible overlap  Total possible overlap  Minimum baseline
...

Coverage: FILENAME_CANDIDATE; timestamps not native-validated.
```

If no eligible triads exist, distinguish:

- `NO_THREE_SITE_CANDIDATE_OVERLAP_WITH_COMPLETE_INVENTORY`
- `THREE_SITE_OVERLAP_BUT_BASELINE_FILTER_EXCLUDES_ALL`
- `OVERLAP_UNDETERMINED_DUE_TO_INCOMPLETE_INVENTORY`
- `NO_ELIGIBLE_ARCHIVE_STATIONS`

Do not say “no telescopes were observing” based only on S3 absence.

## 11. Link back into GNSS discovery pipeline

The network-wide scan should work **without** any `--visibility-results` input. This is essential: Stage-1 results are tied to a particular station triad.

For a selected triad and interval from `network_triad_segments.csv`:

1. Run `tart-common-gnss --stations <id1> <id2> <id3> --start <segment start> --end <segment end> --min-elevation 20 --step 10`, with a contemporaneous GNSS TLE catalogue and correctly handled boundary margins.
2. Run the existing **explicit-stations** `tart-archive-scan --stations ... --visibility-results <new Stage-1 JSON>` to intersect in-band GNSS predictions with candidate archived coverage.
3. Run `tart-retrieve --dry-run` on the selected Stage-1 rank to check exact inventory in the target window.
4. Download native files and invoke Stage-2 validation / `tart2ms` only for viable candidates.
5. Run Stages 3–7 and `tart-campaign` once independent validated observations exist.

Provide **ready-to-copy Stage-1 command templates** for the top triads, but never assert their GNSS satellite visibility until Stage 1 actually runs. Mark templates as proposed next actions, not completed checks.

Do not use an old Stage-1 JSON for a different triad. Do not use current TLEs backward over long historical intervals as if their prediction accuracy were guaranteed.

**Optional follow-on:** Add an export of a campaign plan containing triad IDs, overlap start/end, planned TLE source and proposed Stage-1 invocations; no automated downloads or Stage-3 runs in this feature.

## 12. Geometric nuance and scientific integrity

A terrestrial baseline threshold is a coarse site-diversity filter, **not** the triangulation quality metric. A long baseline can still produce near-parallel sightlines, and GNSS satellite orbital position and station visibility matter.

When ranking scientific candidates after Stage 1, inspect minimum crossing angle, condition number, station elevations, TLE age and in-band signal relevance. Do not label geospatially separated archive candidates “good satellite triangulation” before evaluating these.

A multi-station time interval in the archive may contain **no detectable GNSS signal**. Stage-3 remains the measurement stage.

## 13. Failure handling

- Unknown alias, invalid coordinates, missing archive name: report station exclusion or invalid explicit input; never synthesize station coordinates.
- S3 listing errors: record `ERROR`, retry according to existing robust policy, and preserve uncertainty.
- Partial/stale inventory: allow a useful partial result, but return a non-success quality code and flag it.
- All stations empty or no eligible triads: write empty tables and an explanatory summary; do not crash.
- Missing/invalid cache: invalidate and relist unless `offline`, where unknown days must remain unknown.
- Bad date ranges, negative thresholds, `min_stations` > eligible geographic sites: clear actionable error.
- Fully overlapping co-located stations: never inflate independent-site count.
- Cross-day and sub-second boundaries: deterministic UTC handling.
- No network in unit tests.

Suggested CLI exit behavior aligned with existing scanner:

```text
0 complete scan and report (even when zero candidates)
3 completed but inventory incomplete; scientifically indeterminate absence
2 invalid arguments/configuration
1 unexpected failure
```

## 14. Test specification

Create focused tests in `tests/archive_scan/` and preserve current tests.

### Selection and geometry

- Dynamic station discovery, eligible archive-prefix filtering, alias handling.
- Three geographically distinct stations qualify; three IDs colocated at the same coordinates do not.
- Verified `archive_name` does not imply any files are present.
- 500-km **all-pairs** baseline threshold correctly includes/excludes triads.
- Ellipsoidal geodesic baseline numerical reference test (independent known distances).
- `--all-stations` conflicts with `--stations` and incompatible Stage-1 JSON are handled.

### Inventory and caching

- Count listing calls by unique station/day, *not* by number of triads.
- Full cache reuse, offline cache, expired cache, corrupted/checksum-mismatched cache.
- Day-boundary and prior-day overlap; malformed filenames; duplicate S3 keys.
- Mixed complete, empty, error and unknown station-days.
- Filename last-sample timestamp vs S3 upload timestamp.

### Sweep and interval logic

- Active-station count changes in a synthetic five-station timeline.
- Three stations briefly overlap but for less than 300 s: excluded.
- Triad remains active while a fourth station joins/leaves: **one continuous triad segment**, not fragmented.
- Distinct triads may share the same time interval; network total time must not double count.
- Half-open boundary touch creates zero overlap.
- Genuine file gap is preserved; gap tolerance must not invent covered seconds.
- More than three active sites enumerates only eligible distinct-site triads.
- Deterministic ranking regardless of arbitrary S3 object iteration order.

### Failure and outputs

- One station-day S3 failure does not become a zero-file day.
- No overlap with complete inventory is labelled differently from incomplete coverage.
- No eligible stations produces a clear explanatory result.
- Output CSV/JSON/manifest keys, output SHA-256, metadata sources and hashes.
- Generated Stage-1 command templates use the correct actual triad and interval.
- Existing explicit-station scan behaviour remains byte/semantically compatible where appropriate.
- `tests/test_cli_names.py` still validates all public commands.
- Full unit and regression tests pass without network.

Add a synthetic demo with at least five stations, one colocated pair, time-varying active sets, a station archive outage, and baseline-filtered triads. Verify mathematically expected overlap durations.

## 15. Acceptance criteria

- [ ] `tart-archive-scan --all-stations` implemented and documented.
- [ ] `--min-stations 3` enforces distinct physical sites, not merely three registry IDs.
- [ ] `--min-baseline-km` implements documented minimum all-pair geodesic distance.
- [ ] Station inventory/cache reused; no station/day listed separately for each candidate triad.
- [ ] Output includes station/day availability, station summaries, network active segments, ranked triads and per-triad continuous intervals.
- [ ] Candidate intervals never represented as HDF5-validated observation time.
- [ ] Sweep-line calculation matches synthetic ground truth, including fourth-station join/leave.
- [ ] No fabricated coverage across gaps or day boundaries.
- [ ] Incomplete archive inventory handled explicitly.
- [ ] Original `tart-archive-scan --stations ... --visibility-results ...` still works.
- [ ] Handoff templates for `tart-common-gnss` use new triad-specific Stage-1 runs.
- [ ] Full test suite passes, and counts are reported rather than invented.
- [ ] README/docs/examples, implementation notes, decisions/history log and handoff status updated.

## 16. Repository work management and handoff

Follow the project's engineering standards:

1. Create a dedicated issue/ticket with scope, acceptance tests, and known scientific limitations.
2. Work on a feature branch; record base SHA and final SHA.
3. Record implementation notes and design decisions (event sweep, physical-site grouping, candidate-vs-verified distinction, caching and baseline policy).
4. Keep `docs/handoff/` (or existing handoff convention) current: entry point, commands, tests, assumptions, open issues and follow-up work.
5. Update existing docs rather than creating duplicate contradictory guidance.
6. Document all public Python functions with clear purpose, parameters, return values and algorithm notes; use `logging.info` for key progress and outcomes.
7. Do not commit large observational HDF5 or Measurement Sets.

### Required coding-agent completion report

```text
Branch/base commit/final commit:
Modules changed:
New CLI flags:
Number of eligible archive stations and physical sites:
Inventory station/day prefixes listed / cache hits / errors:
Number of active network segments with >=3 independent sites:
Number of eligible triads before/after baseline filter:
Top triads and their candidate overlap intervals:
Coverage/provenance caveats:
Synthetic tests passed/failed/skipped:
Full tests passed/failed/skipped:
Open issues and next validation commands:
```

**Definition of done:** The user can discover *when and where the archive suggests that three or more geographically distinct TART stations may have recorded simultaneously*, without manually guessing station combinations, while retaining full provenance and honest separation between candidate filename overlap and validated scientific observations.
