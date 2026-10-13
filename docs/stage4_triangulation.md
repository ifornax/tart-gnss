# Stage 4 — multi-station 3-D localisation (`tart-triangulate`)

Stage 4 turns the per-station directions measured by Stage 3
(`tart-localise`) into a 3-D satellite position in ITRS/ECEF, with a full
3×3 covariance, from the geometry of the measured sight lines alone.
TLE predictions never enter the fit: Stage 4 does not read the TLE, and the
Stage-3 TLE columns are ignored. The satellite block of the Stage-3 manifest is
copied to the Stage-4 manifest for later comparison (Stage 6) only.

```text
tart-localise  ->  Az_i(t), El_i(t), cov  (per station, Stage 3)
tart-triangulate -> x(t), y(t), z(t), 3x3 covariance  (ITRS/ECEF, metres)
```

Stage 4 does no orbit fitting, filtering, smoothing, dynamics, TLE residual
analysis or station-bias calibration (Stages 5 and 6).

## 1. Usage

```bash
pip install -e ".[stage3]"          # numpy, pandas, pyarrow, matplotlib
tart-triangulate \
    --stage3-manifest results/stage3/20260625/NORAD_40534/stage3_manifest.json \
    --output-dir results/stage4
```

| Option | Default | Meaning |
|---|---|---|
| `--stage3-manifest` | required | `stage3_manifest.json` written by `tart-localise` |
| `--output-dir` | `results/stage4` | output root; results go to `<root>/<YYYYMMDD>/NORAD_<id>/` |
| `--time-match-tolerance` | 0.75 s | max \|t_station − t_reference\| |
| `--min-stations` | 3 | stations required for a valid estimate (≥ 2) |
| `--method` | `weighted` | `weighted` (Stage-3 covariance) or `unweighted` (closed form) |
| `--bootstrap N` | 0 | parametric bootstrap resamples per epoch |
| `--seed` | 20260625 | bootstrap random seed (recorded) |
| `--stations` | all | restrict to these station ids |
| `--quality-policy` | `strict` | Stage-3 input policy, see §6 |
| `--reference-time` | – | triangulate only the epoch nearest this UTC time |
| `--reference-station` | most rows | station supplying the reference epochs |
| `--max-condition-number` | 1e8 | `POOR_GEOMETRY` threshold |
| `--max-normalised-residual` | 5 | `HIGH_LINE_RESIDUAL` threshold |
| `--max-reduced-chi2` | 5 | `HIGH_CHI2` threshold |
| `--time-mismatch-sigma` | 1 | `TIME_MISMATCH` threshold (see §4) |
| `--no-plots`, `--verbose` | | |

Exit codes: **0** at least one valid 3-D estimate; **3** completed but none
valid; **2** invalid input or configuration; **1** unexpected failure.

## 2. Inputs and frozen station coordinates

Stage 4 reads `stage3_manifest.json` and the direction table next to it
(`directions.parquet`, else `directions.csv`). It needs these columns:
`time_utc, station_id, norad_id, az_obs_deg, el_obs_deg, sigma_az_deg,
sigma_el_deg, cov_az_el, l, m, sigma_l, sigma_m, cov_lm, fit_success,
quality_flags`.

Station coordinates come from the run's own provenance:
`stage3_manifest.stations.<id>.{latitude_deg, longitude_deg, altitude_m}`,
which Stage 3 copied from the Stage-2 manifest. If that Stage-2 manifest is
still reachable and its SHA-256 equals the hash Stage 3 recorded, the
coordinates are cross-checked against it, and any disagreement is an error.
The live station registry is never consulted, so a later registry update
cannot change a past result. The source and the outcome of the cross-check
are recorded per station in `stage4_manifest.json`. Coordinates are converted
to ECEF on the WGS84 ellipsoid. No station is hard-coded.

## 3. Sight lines

Conventions are the same as in Stage 3: local East-North-Up, azimuth from
North (0°) through East (90°), elevation 0° at the horizon and 90° at the
zenith:

    s_ENU = (cos E sin A, cos E cos A, sin E) = (l, m, n),   n = sqrt(1 − l² − m²)

ENU vectors are rotated to ECEF with the transpose of the ECEF→ENU rotation
at the station's geodetic latitude and longitude. The unit vector and its
covariance are built from Stage 3's native `(l, m)` and `cov(l, m)`. These
are finite at the zenith, where `σ_az` is undefined. Az/El and `cov_az_el`
(deg²) are used only if `(l, m)` are missing. The angular covariance of a
sight line is carried as a rank-2 3×3 ECEF matrix tangent to the unit sphere:

    Σ_θ = J Σ_lm Jᵀ,   J = R_ENU→ECEF · ∂s_ENU/∂(l, m) = R · [[1, 0], [0, 1], [−l/n, −m/n]]

Tests check the cardinal directions (e.g. at (0°, 0°) North = +Z,
East = +Y, Up = +X). A non-circular test checks that skyfield's own WGS84
positions and alt/az agree with ours to 10⁻⁶°.

## 4. Time matching and timing diagnostics

Station timestamps are independent: each station dumps on its own clock, and
Stage 3 bins each station from its own first timestamp. Rows are never joined
by row number, and nothing is interpolated.

* The **reference station** is the one with the most valid rows (ties are
  broken alphabetically). `--reference-station` overrides the choice.
* Each of its rows defines a reference epoch `t`. Every other station
  contributes its nearest valid row with `|t_i − t| ≤ tolerance`.
* A row is used in at most one epoch. When two epochs want the same row, the
  closer epoch keeps it.
* The contributions table stores `time_utc` (the reference),
  `station_time_utc` and `station_time_offset_s` for every contributing
  station. `positions.max_time_offset_s` is the spread of the matched
  timestamps.

**Motion during the offset.** The angular rate at each station is estimated
from that station's own measured track (never from a TLE). It is the slope of
a straight-line fit to the unit vectors over a 120 s sliding window; a fit
rather than neighbour differencing keeps measurement noise out of the rate.
The motion is `ω_i |t_i − t|`. `TIME_MISMATCH` is raised, on the station and
the epoch, when this motion exceeds `--time-mismatch-sigma` × the station's
1σ angular error. At that point the timing error is no longer negligible
compared with the measurement error. The flag is a warning, not a blocker.

**Light time.** `ρ_i / c ≈ 66–74 ms` for GNSS is recorded per station
(`light_time_s`) and per epoch (`mean_light_time_s`). It is not corrected. The
satellite moves about 0.25 km in that time, far below the position
uncertainty (§9), so retarded-time modelling is not needed in Stage 4.

## 5. Estimators

Station `i` at `r_i` sees the satellite at `x = r_i + ρ_i s_i` with
`ρ_i > 0`. With `P_i = I − s_i s_iᵀ`, the perpendicular offset of a point from
line `i` is `d_i = P_i (x − r_i)`.

**Unweighted (baseline, `--method unweighted`).** Minimises
`Σ |P_i (x − r_i)|²` by solving `(Σ P_i) x = Σ P_i r_i`. The solution is
closed form, so it needs no initial guess and no TLE.

**Weighted (production, default).** In an orthonormal tangent basis
`E_i = [e1_i e2_i] ⊥ s_i`, the transverse offset is `q_i = E_iᵀ (x − r_i)`.
The basis uses the coordinate axis least aligned with `s_i` as its helper, so
it never degenerates at the poles or the zenith. A direction error `δs_i`
moves the line by `−ρ_i δs_i` at the satellite, so

    Cov(q_i) = ρ_i² C_i,   C_i = E_iᵀ Σ_θ,i E_i   (2×2, rad²)
    χ²(x) = Σ_i q_iᵀ W_i q_i,   W_i = (ρ_i² C_i)⁻¹

The weighted fit iterates to convergence:

1. start from the unweighted solution;
2. compute the ranges `ρ_i = s_iᵀ (x − r_i)`;
3. form `W_i`;
4. solve `(Σ E_i W_i E_iᵀ) x = Σ E_i W_i E_iᵀ r_i`;
5. repeat from step 2 until the position moves by less than 1 mm (at most 50
   iterations).

Ranges are floored at 1 m in absolute value when forming weights. A negative
range is reported, not hidden.

**Numerics.**

* Point estimates use `numpy.linalg.solve`. If the normal matrix has
  condition number above 10¹² or is singular, an SVD least-squares solution
  (`lstsq`) is used instead. Which solver ran is recorded in `solver`, and
  the manifest counts how often the fallback was needed.
* The weighted normal matrix is rescaled to O(1) before solving. `A` is never
  explicitly inverted for the point estimate.
* Non-finite directions, malformed shapes and fewer than two lines raise
  explicit errors.
* A singular or non-positive-definite directional covariance makes the
  weighted fit impossible. The epoch then falls back to the unweighted point
  estimate and is flagged `COVARIANCE_SINGULAR`.

## 6. Stage-3 quality policy

| Stage-3 flag | `strict` (default) | `permissive` |
|---|---|---|
| `FIT_FAILED`, `INVALID_GEOMETRY`, `LOW_BASELINE_COUNT` | rejected | rejected |
| `BOUNDARY_SOLUTION`, `LOW_SNR`, `TLE_UNASSOCIATED` | rejected | kept → epoch `STAGE3_QUALITY_WARNING` |
| `HIGH_FLAGGED_FRACTION`, `HIGH_RESIDUAL`, `IMAGE_DISAGREEMENT`, `ZENITH_AZ_UNCERTAIN` | kept → `INPUT_DIRECTION_WARNING` | same |
| `CALIBRATION_WARNING` (or station calibration warnings) | kept → `CALIBRATION_WARNING` | same |

`strict` therefore uses exactly the rows Stage 3 marked `usable`. Rows without
`fit_success` or without a finite direction are always rejected.

`TLE_UNASSOCIATED` is a post-fit *association* decision taken in Stage 3: the
measured source is more than 5° from the target's prediction. Rejecting such
rows decides which measurements describe the target. It does not put the TLE
into the Stage-4 optimisation; the point is computed from the measured lines
only. `permissive` lets the user see what happens without that association
step.

## 7. Uncertainty

**Linear covariance**, first order in the angular errors:

    Σ_x = N⁻¹                         (weighted,   N = Σ E_i W_i E_iᵀ)
    Σ_x = A⁻¹ (Σ_i ρ_i² E_i C_i E_iᵀ) A⁻¹   (unweighted, A = Σ P_i; sandwich form)

Every valid estimate reports `cov_xx … cov_zz` and `sigma_x/y/z_m`.

The approximation assumes three things:

* Gaussian, unbiased, independent station directions with the Stage-3
  covariance;
* exact station coordinates;
* angular errors much smaller than the crossing angles between the lines.

The covariance is flagged `COVARIANCE_SINGULAR` if it is singular, not
positive definite or non-finite.

**Bootstrap** (`--bootstrap N`, seed recorded). Each station direction is
redrawn from its Stage-3 covariance (a Gaussian in the tangent plane), and
the epoch is triangulated again `N` times. `x/y/z_p16/p50/p84` are reported.
`BOOTSTRAP_UNSTABLE` is raised in any of these cases:

* more than 10 % of resamples fail;
* the bootstrap half-width `(p84 − p16)/2` differs from the linear σ by more
  than a factor of 2;
* the bootstrap median moves from the point estimate by more than
  `max(0.25 σ, 3 × the sampling noise of a median)`. This is the signature of
  the non-linear regime (§8).

## 8. Covariance validation (Monte Carlo)

Set-up:

* exact sight lines from a known position (the three real station sites,
  GNSS-like target at about 20 000 km range, crossing angles 3.0–9.3°);
* `n` realisations perturbed by their angular covariances;
* the empirical scatter compared with the mean predicted covariance
  (`tart_gnss.stage4.uncertainty.monte_carlo_check`;
  `tests/stage4/test_synthetic_uncertainty.py`).

| σ_θ | empirical/predicted σ_x, σ_y, σ_z | along mean LOS | ⟨Mahalanobis²⟩ (expect 3) | mean error / σ |
|---|---|---|---|---|
| 0.05° | 0.97, 0.98, 0.98 | 0.98 | 2.99 | ≤ 0.03 |
| 0.10° | 0.97, 0.97, 0.97 | 0.97 | 2.99 | ≤ 0.02 |
| 0.15° | 0.97, 0.96, 0.97 | 0.97 | 3.00 | ≤ 0.02 |
| 0.5° | 0.98, 0.99, 0.98 | 0.98 | 3.24 | ≤ 0.18 |
| 1.0° | 0.93, 0.94, 0.93 | 0.93 | 3.9 | ≤ 0.31 |
| 2.0° | 0.77, 0.79, 0.76 | 0.77 | 7.0 | ≤ 0.60 |

(n = 1000, weighted estimator; the unweighted sandwich covariance behaves the
same.)

Up to about 0.5° the linear covariance is accurate to a few per cent. From
about 1° the estimator becomes biased: range scales as 1/crossing angle, so
the error distribution is skewed along the line of sight. The predicted σ
also depends on the noisy fitted ranges. This is the regime the bootstrap
median-shift test flags. Stage-3 random errors on real data are 0.06–0.15°,
well inside the linear regime.

The error ellipsoid is long and thin, along the mean line of sight. With
baselines of 1 000–3 000 km against about 20 000 km range, the along-LOS σ
is roughly `ρ σ_θ / sin γ`, about 5–10× the transverse σ `ρ σ_θ`. At
σ_θ = 0.1°, that is about 35 km across the line of sight and 250–300 km
along it.

## 9. Geometry, residual and quality outputs

For each epoch Stage 4 reports:

* **Geometry:** the condition number and rank of `A = Σ P_i`, and the
  minimum and maximum pairwise crossing angle `γ_ij = arccos(s_i·s_j)`.
* **Fit quality:** `chi2`, `dof = 2N − 3` and `reduced_chi2`.
* **Residuals:** per station, the line residual `|P_i(x̂ − r_i)|` (m) and the
  normalised residual `sqrt(q_iᵀ W_i q_i)`; per epoch, the RMS and maximum
  line residual.
* **Ranges and position:** the station ranges, and the geodetic latitude,
  longitude and height of the solution.

| Flag | Meaning | Blocks `valid` |
|---|---|---|
| `INSUFFICIENT_STATIONS` | fewer than `--min-stations` matched; 2-station solutions are written as degraded diagnostics | yes |
| `TIME_MISMATCH` | target motion during a station's time offset > `--time-mismatch-sigma` × σ_θ | no |
| `POOR_GEOMETRY` | condition number > `--max-condition-number` (1e8) or rank < 3 | yes |
| `NEGATIVE_RANGE` | some `ρ_i ≤ 0` (point behind a station) | yes |
| `HIGH_LINE_RESIDUAL` | a station's normalised residual > 5 | no |
| `HIGH_CHI2` | reduced χ² > 5 | no |
| `INPUT_DIRECTION_WARNING` | a contributing Stage-3 row has an advisory flag | no |
| `CALIBRATION_WARNING` | a contributing station has a calibration warning | no |
| `STAGE3_QUALITY_WARNING` | a Stage-3 blocking flag let through by `permissive` | no |
| `COVARIANCE_SINGULAR` | directional or position covariance unusable | yes |
| `BOOTSTRAP_UNSTABLE` | see §7 | no |
| `TRIANGULATION_FAILED` | no finite, full-rank, converged solution | yes |

A **valid** epoch is the input for Stage 5. Warnings mark epochs that need
down-weighting or inspection.

## 10. Outputs

```text
results/stage4/<YYYYMMDD>/NORAD_<id>/
├── positions.{parquet,csv}             one row per reference epoch
├── station_contributions.{parquet,csv} one row per station per epoch
├── stage4_manifest.json
├── diagnostics/  position_track_ecef.png  geometry_condition.png  station_residuals.png
│                 position_uncertainty.png timing_offsets.png
└── logs/stage4.log
```

**`positions`** has:

* the required §21 fields: `time_utc, time_mjd, norad_id, x/y/z_ecef_m,
  cov_xx…cov_zz, sigma_x/y/z_m, n_stations, rms/max_line_residual_m,
  condition_number, min/max_crossing_angle_deg, chi2, reduced_chi2,
  fit_success, quality_flags, latitude_deg, longitude_deg, altitude_m`;
* `reference_station, method, stations, matrix_rank, dof,
  max_normalised_residual`;
* `max_time_offset_s, max_motion_during_offset_deg, mean_light_time_s`;
* the bootstrap percentiles;
* `solver, n_iterations, converged, valid`.

**`station_contributions`** has:

* the required §22 fields;
* `l, m, sigma_theta_deg, direction_source`, and the ECEF sight vector;
* `normalised_residual`;
* `angular_rate_deg_s, motion_during_offset_deg, light_time_s`.

Flags are lists in Parquet and `;`-joined in CSV.

**`stage4_manifest.json`** (schema 1) records:

* the input manifest and direction table with their SHA-256 hashes;
* the frame and units;
* `time_matching`: tolerance, reference station rule, offsets, rate source,
  light-time policy;
* `triangulation`: method, estimators, solver and fallback count,
  `tle_used: false`;
* `uncertainty`: model, bootstrap count, seed;
* `quality_policy` with rejection counts, `quality_thresholds` and the
  blocking flags;
* a summary with flag counts;
* per-station coordinates, source, cross-check and residual statistics;
* the satellite reference block (copied, unused);
* software versions and git state at run start;
* provenance: command line, time, working directory, the Stage-3 command;
* the output paths.

## 11. Synthetic end-to-end demonstration

`python scripts/stage4_synthetic_demo.py` writes
[`docs/examples/stage4`](examples/stage4).

* **Truth.** One 3-D track: the TLE ITRS track of PRN 26 on 2026-06-25, moved
  by a fixed **1000 km perpendicular to the mean line of sight**.
* **Simulation.** Three synthetic 24-antenna stations at the real station
  coordinates observe the truth, with timestamp offsets of 0, +0.25 and
  −0.40 s. Complex visibilities are simulated, `tart-localise` measures each
  station's Az/El (1 s bins), and `tart-triangulate` intersects the lines.

Middle epoch (from `demo_report.txt`):

| | x (km) | y (km) | z (km) |
|---|---|---|---|
| true | 16151.034 | 7563.587 | −19951.461 |
| recovered | 16216.264 | 7604.119 | −19985.951 |
| error | 65.2 | 40.5 | −34.5 |
| predicted σ | 228.7 | 88.7 | 324.4 |

* **3-D error:** 84 km.
* **Distance from the TLE position:** 989 km.
* **Geometry:** condition number 261; crossing angles 2.59–8.20°; reduced χ²
  0.95; time-offset spread 0.625 s.

| Station | offset (s) | range (km) | line residual (km) | normalised | σ_θ (°) |
|---|---|---|---|---|---|
| madagascar (reference) | 0.000 | 21884 | 54.1 | 1.25 | 0.143 |
| mauritius | −0.250 | 22194 | 48.3 | 0.90 | 0.172 |
| south_africa | +0.375 | 20749 | 35.4 | 0.68 | 0.148 |

All 30 epochs:

* **Validity:** 29 valid. The last has only two stations, because of the
  offsets.
* **3-D error:** median 293 km, dominated by the along-line-of-sight
  component (median 290 km); transverse median 32 km.
* **Pulls (x, y, z):** std 0.95, 0.87, 0.99.
* **Recovered displacement from the TLE** along the injected offset
  direction: median **996 km** (injected 1000 km; scatter 35 km).

The position is therefore derived from the measured directions, not from the
TLE. Stage 4 never reads the TLE. A test (`test_tle_never_used`) shows the
result is bit-identical with the TLE columns and satellite block removed.

## 12. Real-data run (2026-06-25, GPS PRN 26, NORAD 40534)

Stage 3 was re-run with the current code on the real three-station Stage-2
data (`tart-localise --all-stations --image-validation`). Stage 4 was then run
on its output. The example lives in
[`docs/examples/stage4_real`](examples/stage4_real), with its 10 s Stage-3
input. Large real datasets (Measurement Sets, HDF5) are not committed.

**Spec acceptance command, default options.**

* **1 s Stage-3 bins** (the Stage-3 default): exit 3. Only 1 of 1980 Stage-3
  rows is usable (per-bin SNR 3–4). Strict rejection: 1166 `LOW_SNR`, 813
  `BOUNDARY_SOLUTION`. No epoch has three stations.
* **10 s Stage-3 bins:** exit 3.
  * 107 of 204 rows are usable. Strict rejection: 86 `BOUNDARY_SOLUTION`,
    11 `TLE_UNASSOCIATED`. ZA-Rhodes has 68 usable rows, MU-UDM 37 (from
    09:18:36), MG-Tana 2.
  * Each station's bins start at its own first timestamp, so MU's bin centres
    sit 2.6 s from ZA's in every bin. The 0.75 s tolerance never pairs them,
    and no epoch has three stations.

**Diagnostic run** (`--time-match-tolerance 6 --min-stations 2 --bootstrap 200`;
the tolerance is widened and two-station epochs are accepted, both
deliberately and for diagnosis only):

| | value |
|---|---|
| matched epochs with ≥ 2 stations | 36 (35 ZA+MU, 09:18:34–09:24:22; 1 ZA+MU+MG at 09:23:53.9) |
| time offsets | MU 2.61–2.62 s; MG 0.26 s; target motion during the offset 0.01–0.05° (< σ_θ, so no `TIME_MISMATCH`) |
| crossing angles / condition number | ZA–MU 7.4–7.7°; MG–MU 3.9° · 222–306 |
| Stage-3 σ_θ | 0.14–0.16° |
| position σ (x, y, z) | median 264, 165, 456 km; bootstrap agrees (no `BOOTSTRAP_UNSTABLE`) |
| line residuals | median 116 km (ZA), 125 km (MU); MG 448 km in the 3-station epoch |
| reduced χ² | median 11.5 (26 of 36 epochs `HIGH_CHI2`) |
| solution | geodetic height 22 690–23 880 km (geocentric radius ≈ 29 700 km) |
| continuity | median step 130 km between 10 s epochs (σ-dominated; the satellite moves about 39 km) |

The single three-station epoch (09:23:53.9) gives x, y, z =
15 499, 10 245, −23 315 km, with σ = 228, 146, 414 km. The flags are
`HIGH_LINE_RESIDUAL` (MG, normalised residual 7.4) and `HIGH_CHI2`
(reduced χ² 24).

**Interpretation, without TLE comparison (Stage 6).**

* The solutions sit about 3 000 km further out along the line of sight than
  a GPS orbit (about 20 200 km height) would allow. The reduced χ² is much
  greater than 1.
* Both are what the Stage-3 systematics predict. The directions carry slowly
  varying offsets of 1–2°, much larger than their 0.15° random σ. With a
  7° crossing angle, a 1° relative direction error moves the intersection by
  about 21 000 km × 1°/7°, i.e. about 3 000 km along the line of sight.
* The geometry is good (condition number about 240, rank 3, no
  `POOR_GEOMETRY`). The limitation is the accuracy of the Stage-3 directions
  and the bin alignment, not the triangulation.
* Stage 4 does not calibrate station biases (by design).

**What would make the default run produce three-station estimates:**

* common UTC-aligned Stage-3 bin edges, so MU and ZA bins coincide;
* more usable MG rows, since MG-Tana is dominated by a non-target source;
* calibration of the Stage-3 direction systematics (Stage 5 or later).

## 13. Limitations

* **Stage-3 systematics are not modelled.** On real data the slowly varying
  Stage-3 direction systematics (0.1–0.2° wander, 1–2° offsets) exceed the
  random σ. Positions will show `HIGH_CHI2` and residuals that persist over
  minutes. Stage 4 reports this rather than inflating σ: per-station bias
  calibration is deferred by design.
* **Epochs are independent.** Successive epochs share those systematics, so
  they are not independent measurements. Stage 5 should allow for correlated
  errors.
* **No light-time or aberration correction.** Both are below the
  uncertainty, as §4 explains.
* **Long, thin uncertainty.** The geometry (baselines much smaller than the
  range) makes the uncertainty strongly elongated along the line of sight,
  whatever the estimator.
* **Bin alignment.** Stage 3 starts each station's bins at that station's
  first timestamp. With bins wider than the tolerance, stations can therefore
  fail to match at all (§12). Aligning Stage-3 bins to a common UTC grid
  would fix this; it is a Stage-3 change and is left as a follow-up.

## 14. CLI naming convention

Every public command of this project starts with `tart-`. The canonical
commands are:

    tart-stations-update  tart-stations-list  tart-common-gnss  tart-retrieve  tart-localise  tart-triangulate

Future stages must follow the same pattern. `tests/test_cli_names.py` checks
that:

* the entry points in `pyproject.toml` are exactly these;
* each `argparse` `prog=` matches its entry-point name;
* no legacy name (`common-gnss`, `retrieve-tart`, `update-tart-stations`,
  `list-tart-stations`) appears in active code or docs.

Historical provenance files that record what was actually executed are never
rewritten.
