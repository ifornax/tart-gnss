# Implementation Specification — Dynamic TART Station Registry Integration

## 1. Objective

Extend the existing `tart-gnss` pipeline so that TART station metadata can be retrieved dynamically from the official TART network registry:

```text
https://map.elec.ac.nz/api/v1/telescopes
```

The current repository already has:

- `src/tart_gnss/stations.py`
- `src/tart_gnss/config/stations.yaml`
- Stage-1 CLI: `common-gnss`
- Stage-2 CLI: `tart-retrieve`
- packaged station metadata containing `tart_id`, `api_url`, `archive_name`, coordinates, and provenance.

The new implementation MUST integrate with these existing abstractions rather than creating a parallel station system.

The main goal is to allow the pipeline to:

1. discover all TART deployments currently published by the registry;
2. validate and normalize their coordinates;
3. preserve metadata provenance;
4. generate/update a machine-readable station registry;
5. let Stage 1 and Stage 2 use any valid TART station without manual code changes;
6. preserve reproducibility by distinguishing live registry metadata from frozen run metadata.

---

## 2. Current Repository State

The current `Station` dataclass already contains:

```python
id
name
latitude_deg
longitude_deg
altitude_m
source
tart_id
api_url
archive_name
```

Do NOT replace this model unless strictly necessary.

The current static configuration:

```text
src/tart_gnss/config/stations.yaml
```

already records manually curated aliases and authoritative coordinates. The implementation should evolve the existing station layer from:

```text
manually maintained YAML
```

to:

```text
live registry -> normalized registry snapshot -> validated Station objects
```

while retaining support for manually curated aliases.

---

## 3. Authoritative Registry Source

Use:

```text
https://map.elec.ac.nz/api/v1/telescopes
```

as the primary network-discovery endpoint.

The implementation MUST NOT scrape the interactive map HTML.

Treat the registry response as external data whose schema may evolve. Therefore:

- isolate parsing in one module;
- validate required fields explicitly;
- fail clearly if the schema changes;
- preserve the raw downloaded payload;
- record retrieval time;
- record a SHA-256 checksum.

Do not silently guess missing coordinates.

---

## 4. Design Principle

Separate these concepts:

### 4.1 Registry identifier

The official TART telescope identifier, e.g.

```text
mu-udm
mg-tana
za-rhodes
```

This should remain the canonical machine identifier.

### 4.2 User-friendly alias

Examples:

```text
mauritius
madagascar
south_africa
```

Aliases are convenience labels and may remain manually curated.

### 4.3 Display name

Human-readable station/site name from TART metadata.

### 4.4 Registry snapshot

The exact station metadata downloaded at a known UTC time.

Do not use aliases as the authoritative network identity.

---

## 5. Recommended Architecture

Add:

```text
src/tart_gnss/
├── registry.py
├── stations.py
└── config/
    ├── stations.yaml
    └── station_aliases.yaml
```

Optionally add:

```text
data/stations/
```

for timestamped registry snapshots.

Responsibilities:

### `registry.py`

- HTTP retrieval;
- registry JSON parsing;
- schema validation;
- normalization;
- SHA-256;
- provenance;
- writing snapshots;
- comparing snapshots.

### `stations.py`

Remain responsible for:

- `Station` dataclass;
- coordinate validation;
- loading station configurations;
- resolving station identifiers and aliases.

Avoid putting HTTP/network logic into `stations.py`.

---

## 6. Registry Data Model

Introduce a normalized registry record.

```python
@dataclass(frozen=True)
class RegistryStation:
    tart_id: str
    name: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    api_url: Optional[str]
    archive_name: Optional[str]
    registry_source: str
    registry_retrieved_utc: str
```

Keep Python >=3.9 compatibility.

Every public class/function must have a clear docstring describing:

- purpose;
- inputs;
- outputs;
- units;
- failure behaviour.

---

## 7. Registry Fetch Function

Implement something conceptually like:

```python
def fetch_registry(
    url: str = DEFAULT_REGISTRY_URL,
    timeout_s: float = 30.0,
) -> dict:
    """Download the current TART telescope registry.

    Parameters
    ----------
    url
        TART network registry endpoint.
    timeout_s
        HTTP request timeout in seconds.

    Returns
    -------
    dict
        Parsed registry JSON.

    Raises
    ------
    RegistryFetchError
        If the endpoint cannot be reached or returns invalid JSON.
    """
```

Use a real HTTP client such as `requests` or standard-library `urllib`.

Prefer adding:

```text
requests>=2.28
```

to the appropriate dependency group if not already present.

Use `logging.info` for fetch status.

---

## 8. Preserve Raw Registry Snapshot

Every successful refresh should preserve the exact raw registry payload.

Recommended layout:

```text
data/stations/
├── registry_20261009T154800Z.json
├── registry_20261009T154800Z.meta.json
└── latest.json
```

Metadata sidecar:

```json
{
  "source_url": "https://map.elec.ac.nz/api/v1/telescopes",
  "retrieved_utc": "...",
  "sha256": "...",
  "station_count": 0
}
```

Do not overwrite historical snapshots.

---

## 9. Registry Parser

Implement a single parser such as:

```python
def parse_registry(payload, retrieved_utc):
    """Normalize the TART registry into validated station records."""
```

The parser must:

1. iterate over every published telescope;
2. extract the canonical TART ID;
3. extract coordinates;
4. validate coordinates using the existing validation logic;
5. record malformed entries separately;
6. return valid station records plus diagnostics.

Do not scatter nested registry JSON access throughout the codebase.

Do not fail the entire registry merely because one station is malformed unless strict mode is enabled.

---

## 10. Coordinate Validation

Reuse the existing:

```python
validate_coordinates(...)
```

from `stations.py`.

Do NOT duplicate latitude/longitude/altitude validation rules.

Reject or flag:

- missing coordinates;
- NaN/infinite values;
- latitude outside `[-90, 90]`;
- longitude outside `[-180, 180]`;
- implausible altitude outside the existing accepted range.

Do not substitute `0` for missing values.

---

## 11. API URL Generation

If the registry provides an API URL, use it.

If it provides only a TART ID and the current project intentionally uses the TART API convention, generate:

```text
https://api.elec.ac.nz/tart/<tart_id>
```

only in a dedicated function:

```python
def default_api_url(tart_id: str) -> str:
    ...
```

Document that this is convention-derived metadata, not registry-originated metadata.

---

## 12. Archive Name Handling

Do NOT assume every registry telescope has archived observation data.

Keep:

```text
registry station existence != archive existence
```

Recommended fields:

```text
archive_name
archive_verified
archive_verified_utc
```

If archive discovery verifies a prefix, record it.

If archive availability has not been checked:

```text
archive_name: null
archive_verified: false
```

Do not mark a station archive-capable merely because it exists in the registry.

---

## 13. Alias Configuration

Move convenience aliases into a manually curated mapping.

Example:

```yaml
aliases:
  mauritius: mu-udm
  madagascar: mg-tana
  south_africa: za-rhodes
  south_africa_hammanskraal: za-hammanskraal
  namibia: na-unam
  zambia: zm-cbu
  ghana: ghana
```

Canonical identity remains the TART ID.

Both should work:

```bash
common-gnss --stations mu-udm mg-tana za-rhodes ...
```

and:

```bash
common-gnss --stations mauritius madagascar south_africa ...
```

---

## 14. Station Resolution Changes

Update `resolve_stations()` so it resolves in this order:

1. exact canonical station ID;
2. normalized canonical station ID;
3. configured alias;
4. `tart_id` where relevant.

Continue to reject duplicates.

Examples that should resolve to the same station:

```text
za-rhodes
za_rhodes
south_africa
```

Do not silently break historical Stage-1 output schemas. If canonical IDs change in output, version the schema.

---

## 15. Registry-to-Station Conversion

Provide a function such as:

```python
def registry_to_stations(
    registry_stations,
    aliases=None,
    archive_metadata=None,
):
    """Convert normalized registry entries into validated Station objects."""
```

Preserve:

- `tart_id`;
- coordinates;
- display name;
- API metadata;
- archive metadata;
- source/provenance.

Resulting objects must remain compatible with Stage 1 and Stage 2.

---

## 16. CLI: Refresh Registry

Add a new CLI entry point:

```text
tart-stations-update
```

Example:

```bash
tart-stations-update
```

Optional arguments:

```text
--registry-url URL
--output PATH
--snapshot-dir PATH
--aliases PATH
--verify-archives
--strict
--dry-run
--verbose
```

Expected summary:

```text
Fetching TART registry...
Found N registry entries.
Validated M stations.
K entries skipped due to invalid/missing coordinates.
Wrote registry snapshot ...
Wrote normalized station catalogue ...
```

Never fabricate a station to make counts match expectations.

---

## 17. CLI: List Stations

Add:

```text
tart-stations-list
```

Example output:

```text
ID              Alias          Site                    Lat        Lon       Alt(m)
--------------------------------------------------------------------------------
mu-udm          mauritius      ...                     ...        ...       ...
mg-tana         madagascar     ...                     ...        ...       ...
za-rhodes       south_africa   ...                     ...        ...       ...
```

Useful options:

```text
--archive-only
--format table
--format json
--format yaml
```

Only add `--online-only` if online state can be verified reliably.

---

## 18. Integrate With `common-gnss`

Stage 1 should support frozen and refreshed catalogues.

Recommended option:

```text
--refresh-stations
```

However, for reproducibility, live refresh MUST NOT become the hidden default for every Stage-1 run.

Preferred behaviour:

```text
default -> use packaged/frozen normalized registry
--refresh-stations -> explicitly fetch latest registry first
```

---

## 19. Run Provenance

Stage-1 metadata must record the exact station catalogue used.

Add:

```text
station_registry_source
station_registry_retrieved_utc
station_registry_sha256
station_catalogue_path
```

The three selected stations must still be copied into run metadata with exact coordinates.

A future registry update must not change the interpretation of an old run.

---

## 20. Integrate With Stage 2

Stage 2 should use the exact station coordinates embedded in Stage-1 output wherever possible.

Do NOT silently refresh station coordinates during Stage 2.

Stage 2 may consult the current registry for operational metadata, but if current coordinates differ from the frozen Stage-1 values, report:

```text
station_metadata_changed: true
```

and retain both versions.

---

## 21. Registry Diff Capability

Implement a diff capability:

```bash
tart-stations-update --diff
```

or equivalent.

Detect:

```text
new stations
removed stations
coordinate changes
name changes
metadata changes
```

For coordinate changes, report numerical deltas.

Do not silently overwrite scientifically used station metadata without recording the change.

---

## 22. Optional Geographic Distance Check

For coordinate updates, calculate approximate physical displacement using an established geodesy/Astropy transformation.

Report displacement in metres.

Large changes should trigger a warning.

Do not invent a universal rejection threshold without justification.

---

## 23. Offline Behaviour

Support offline reproducibility.

If live retrieval fails:

- do not destroy existing station data;
- do not silently substitute another source;
- clearly report which frozen snapshot remains available.

Recommended:

```text
tart-stations-update
```

fails if the live registry cannot be retrieved.

But:

```text
common-gnss
```

continues using its existing frozen station catalogue unless `--refresh-stations` was explicitly requested.

---

## 24. Error Types

Add explicit errors:

```python
class RegistryError(Exception):
    """Base error for TART registry operations."""

class RegistryFetchError(RegistryError):
    """Raised when the remote TART registry cannot be retrieved."""

class RegistrySchemaError(RegistryError):
    """Raised when the registry payload has an unsupported schema."""

class RegistryStationError(RegistryError):
    """Raised when an individual registry station is invalid."""
```

Do not use bare `except:`.

---

## 25. Logging

Use `logging`.

Examples:

```python
logging.info("Fetching TART registry from %s", registry_url)
logging.info("Registry contains %d telescope entries", count)
logging.warning("Skipping station %s: missing valid coordinates", tart_id)
```

Operational diagnostics should not use `print()`.

---

## 26. Dependencies

If using `requests`, update `pyproject.toml`, e.g.

```toml
dependencies = [
    ...
    "requests>=2.28",
]
```

Keep Python >=3.9 compatibility.

Do not introduce a heavy framework for a simple HTTP JSON endpoint.

---

## 27. Recommended Files

Add:

```text
src/tart_gnss/registry.py
src/tart_gnss/config/station_aliases.yaml
tests/test_registry.py
```

Optional wrapper:

```text
scripts/refresh_station_registry.py
```

Potential generated packaged file:

```text
src/tart_gnss/config/stations_registry.yaml
```

Clearly distinguish generated station data from hand-maintained aliases.

---

## 28. Generated Normalized Catalogue

Example schema:

```yaml
metadata:
  source: "https://map.elec.ac.nz/api/v1/telescopes"
  retrieved_utc: "..."
  sha256: "..."

stations:

  mu-udm:
    name: "..."
    tart_id: "mu-udm"
    latitude_deg: ...
    longitude_deg: ...
    altitude_m: ...
    api_url: "https://api.elec.ac.nz/tart/mu-udm"
    archive_name: null
    source: "TART network registry snapshot ..."

  mg-tana:
    ...
```

Do not replace unknown fields with guesses.

---

## 29. Tests

Tests are mandatory and should avoid live-network dependence.

Create local registry fixtures.

Required tests:

### Valid registry
All valid entries become normalized stations.

### Missing coordinates
`null` latitude must never become `0`.

### Invalid longitude
Longitude > 180 must be rejected by existing coordinate validation.

### Duplicate TART ID
Duplicate canonical IDs must raise or be explicitly diagnosed.

### Alias resolution
Verify `mauritius` and `mu-udm` resolve to the same station.

### Unknown alias
Fail clearly.

### Snapshot hashing
Verify deterministic SHA-256.

### Registry diff
Test station added, removed, coordinate changed, and name changed.

### Fetch failure
Mock timeout/HTTP failure and confirm existing local catalogue is untouched.

### Schema change
Unexpected top-level structure must raise `RegistrySchemaError`.

---

## 30. Backward Compatibility

Existing commands must continue to work:

```bash
common-gnss   --stations mauritius madagascar south_africa   ...
```

Do not break Stage 1 or Stage 2 interfaces.

Existing `stations.yaml` semantics should either remain supported or be migrated through a documented compatibility loader.

---

## 31. Scientific Reproducibility Rule

The live registry is for **discovery**.

It must not make old analyses non-reproducible.

Every Stage-1 run must retain the exact:

```text
station coordinates
station IDs
registry retrieval timestamp
registry snapshot hash
```

used for that run.

Future registry refreshes must not alter existing results.

---

## 32. Suggested Workflow

Refresh:

```bash
tart-stations-update
```

List:

```bash
tart-stations-list
```

Choose any three:

```bash
common-gnss   --stations mu-udm mg-tana za-rhodes   --start 2026-10-10T00:00:00Z   --end 2026-10-11T00:00:00Z
```

Aliases remain valid:

```bash
common-gnss   --stations mauritius madagascar south_africa   --start ...   --end ...
```

The Stage-1 result must freeze the selected station metadata.

---

## 33. Acceptance Criteria

Complete when:

1. live registry can be downloaded;
2. raw payload is preserved with timestamp and SHA-256;
3. valid deployments normalize into the existing `Station` model;
4. malformed stations are clearly diagnosed;
5. coordinate validation reuses existing rules;
6. canonical TART IDs are retained;
7. aliases still work;
8. `tart-stations-list` displays the discovered network;
9. `tart-stations-update` creates a reproducible normalized snapshot;
10. Stage 1 can use any three discovered stations;
11. Stage 2 remains compatible with the same station objects;
12. archive availability is not conflated with registry presence;
13. old Stage-1 results remain reproducible;
14. registry snapshots can be diffed;
15. current `mauritius madagascar south_africa` usage remains valid;
16. tests pass without requiring live network access.

---

## 34. Do Not Do

Do NOT:

- scrape the TART map HTML;
- hard-code a new list of coordinates in Python;
- assume every registered TART has archive data;
- assume every registered TART API is online;
- replace missing coordinates with zero;
- silently update coordinates used by prior analyses;
- remove current alias support;
- make live network access mandatory for normal Stage-1 operation;
- mix station-registry retrieval with TLE propagation;
- mix station-registry retrieval with Stage-2 HDF5 retrieval.

The station registry is an infrastructure/provenance component and should remain separate from the scientific visibility and localisation calculations.
