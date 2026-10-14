# Stage 6: TLE trajectory validation and residual analysis (`tart-validate`)

Stage 6 is the first stage where a TLE enters the scientific comparison. It
takes the Stage-5 trajectory, which never saw a TLE: the smoothed GCRS
position and velocity, with the full 6×6 covariance at every epoch. It then
propagates *one selected* element set with SGP4 to the exact Stage-5 epochs
and reports residuals.

**Residual sign convention, permanent:** `observed − TLE`, for every residual
in every table. The tests check it.

The TLE is a *reference*, not ground truth. It has no formal covariance, so
there are two kinds of statistic:

- *measurement-conditional* ones use only the Stage-5 covariance and hold the
  TLE fixed;
- *reference-aware* ones add an **empirical reference dispersion**, taken
  from an ensemble of nearby TLEs. This is never called an
  orbit-determination covariance.

Stage 6 never fits, shifts or adjusts the TLE. It never subtracts a time
offset or an RTN bias, never rescales the Stage-5 covariance, and never
writes to the Stage-5 run. The last point is checked: the Stage-5 manifest
and trajectory are re-hashed after the run, and the result is recorded as
`stage5_inputs_unchanged`. A test hashes every file of the Stage-5 directory
before and after. Attribution
of a residual belongs to Stage 7.

Specification: [`tart_gnss_stage6_implementation.md`](../tart_gnss_stage6_implementation.md).

## 1. Inputs

`stage5_manifest.json` must have `stage == 5`, `schema_version == 1` and
`state_frame == "GCRS"`. Stage 6 also reads `trajectory.{parquet,csv}`,
looking next to the manifest first. Both are hashed.

The 6×6 smoothed covariance is rebuilt in full from its 21 upper-triangle
columns, using `tart_gnss.stage5.output.state_covariance`.

**Stage-5 input policy (`--quality-policy`):**

| | Rows compared |
|---|---|
| `strict` (default) | Stage-5 `valid` rows with a finite smoothed state: none of `NOT_INITIALISED`, `SMOOTHER_FAILED`, `COVARIANCE_INVALID`, `INSUFFICIENT_MEASUREMENTS` |
| `permissive` | any row with a finite state. Where the smoothed state is missing, the filtered state and covariance are used (`stage5_state_source = filtered`) |

Rows that are not compared stay in the table with NaN residuals, the flag
`INVALID_STAGE5_STATE`, and a `stage5_rejection_reason`.

Advisory Stage-5 flags are copied verbatim into `stage5_quality_flags`, and
the row gets `STAGE5_QUALITY_WARNING`. These are `FILTER_DIVERGENCE`,
`INNOVATION_OUTLIER`, `LARGE_TIME_GAP`, `EXTRAPOLATED`, `HIGH_NIS`,
`DYNAMICS_MISMATCH`, `STAGE4_QUALITY_WARNING` and `INITIAL_VELOCITY_WEAK`.

## 2. Reference TLE

Exactly one primary source is required. The two options are mutually
exclusive, and exit 2 if neither or both are given.

| Option | Meaning |
|---|---|
| `--tle-file PATH` | Every set of the target NORAD ID in the file is read; a history file may hold several epochs. If several match, the one whose epoch is nearest the trajectory midpoint is used, and the rule and candidate epochs are recorded. A `.meta.json` sidecar, if present, must carry the file's SHA-256, or the file is rejected as modified. |
| `--tle-source stage1` | The exact Stage-1 element set, recovered through the provenance: Stage-5 manifest → `input_stage4_manifest` → `satellite_reference.satellite`, which Stages 1–3 copy unchanged. The chain and the Stage-4 manifest hash are recorded. |
| `--tle-source download` | An explicit CelesTrak download, archived in `--tle-download-dir` with its provenance sidecar before use. CelesTrak serves only *current* elements, so this is never suitable for past observations; the age checks flag it. Nothing is ever downloaded implicitly. |

Two cases are errors, and both exit 2:

- **`TLE_NOT_FOUND`:** the file has no set for this NORAD ID.
- **`TLE_NORAD_MISMATCH`:** the file holds a single set, for another object.

The manifest records the exact lines, epoch, file, SHA-256, source, selection
rule and the SGP4 implementation (`sgp4` version, via Skyfield, WGS72).

## 3. TLE age

Age is `t_obs − t_TLE epoch`, in days:

- **positive:** forward propagation;
- **negative:** propagation backwards from the epoch (`TLE_BACKPROPAGATED`).

`TLE_AGE_WARNING` is raised above `--tle-age-warning-days` (3), and
`TLE_STALE` above `--tle-age-stale-days` (14). These are pipeline quality
thresholds, not SGP4 accuracy limits. Each row carries
`tle_epoch_utc`, `tle_age_days` and `tle_propagation_direction`.

## 4. Propagation and frames

The TLE is propagated directly with SGP4 at **every** Stage-5 epoch. Nothing
is interpolated.

| Item | Implementation |
|---|---|
| SGP4 source frame | TEME (true equator, mean equinox of date), `sgp4` library, WGS72 |
| transformation library | Skyfield `EarthSatellite.at(t)`: TEME → ITRS with GMST (IAU 1982), then ITRS → GCRS through the inverse of `itrs.rotation_at(t)` (IAU 2000A/2006, GAST from UT1). Velocity is carried through, including the Earth-rotation term, giving the **inertial** GCRS velocity. |
| comparison frame | GCRS, the same rotation Stage 5 uses. ECEF positions of both trajectories are `R(t) r_GCRS` with one shared matrix. |
| time scale | UTC from the Stage-5 table through `tart_gnss.stage5.frames.skyfield_times` (whole days plus seconds of day, builtin timescale, leap seconds handled) |
| Earth orientation | polar motion zero (≲ 40 m at GNSS range); UT1 from Skyfield's bundled tables; no network |

TEME is never treated as GCRS/J2000. In 2026 the two differ by about 0.36°,
or about 170 km at GNSS radius, and a test asserts the difference.

Earth-fixed velocities are never compared with inertial ones. They differ by
|ω×r| ≈ 1.9 km/s, and a test asserts that difference too.

**Frame validation** (`tests/stage6/test_frames.py`) uses an *independent*
route: raw `sgp4` TEME, rotated with an explicit IAU-1982 GMST (Vallado,
eq. 3-47), then the Stage-5 ITRS → GCRS matrices. Skyfield's `EarthSatellite`
is not involved. The two routes agree to **2.5 cm in position and 5 µm/s in
velocity**. Two further checks:

- Skyfield's ECEF position equals `R(t) r_GCRS` to 1 mm.
- The GCRS velocity matches a central difference of the position to a few
  cm/s. SGP4's analytic velocity is not the exact derivative of its
  position.

## 5. RTN frame

The RTN axes (= RIC here) come from the **TLE-predicted** state:

```text
R = r_TLE / |r_TLE|                                radial
N = (r_TLE × v_TLE) / |r_TLE × v_TLE|              cross-track (normal, "C")
T = N × R                                          along-track / in-track ("I")
Q = [Rᵀ; Tᵀ; Nᵀ]      det Q = +1,  R × T = N
```

The tests check orthonormality, handedness, the identity basis of a
canonical circular orbit, and that a degenerate geometry raises
`RTNGeometryError`.

The residuals are expressed in these axes:

- **Position:** `dr_RTN = Q Δr`.
- **Velocity:** `dv_RTN = Q Δv`, the *inertial* velocity difference projected
  on the same axes. It is not the time derivative of the RTN position
  residual, which also contains the frame rotation ω×Δr: about 0.15 m/s per
  km at GNSS orbital rate.
- **Covariance:** the Stage-5 position and velocity blocks are rotated as
  `C_RTN = Q C Qᵀ`, from the full covariance. This gives `sigma_*_obs`,
  `cov_radial_along/radial_cross/along_cross_obs_m2` and `sigma_dv_*_obs`.
  Test 34.5 checks it against rotated Monte Carlo samples and against a
  covariance built in RTN.

**Time-shift signature (test).** Take observed = SGP4(t + τ). Then:

- `dr_T ≈ |v|τ`;
- `dv_R ≈ −GM/r²·τ`, from gravity;
- `equivalent_along_track_time_offset_s = dr_T/|v_TLE| ≈ τ`.

The time offset is a diagnostic. It is never read automatically as a clock
error.

## 6. Residuals and distances (per epoch)

| Column | Definition |
|---|---|
| `dr_{x,y,z}_m`, `position_residual_3d_m` | Δr = r_obs − r_TLE (GCRS), \|Δr\| |
| `dv_{x,y,z}_m_s`, `velocity_residual_3d_m_s` | Δv = v_obs − v_TLE (inertial GCRS), \|Δv\| |
| `dr_{radial,along_track,cross_track}_m`, `dv_*_m_s` | RTN projections (§5) |
| `geocentric_angular_residual_deg` | angle between r_obs and r_TLE (atan2 form, accurate at small angles) |
| `geocentric_range_residual_m` | \|r_obs\| − \|r_TLE\| |
| `equivalent_along_track_time_offset_s` | Δr_T / \|v_TLE\| |
| `position_mahalanobis_obs_cov` | D² = Δrᵀ C_r⁻¹ Δr with the Stage-5 covariance only, **conditional on the TLE being fixed** |
| `state_mahalanobis_obs_cov` | D² = Δxᵀ P⁻¹ Δx (6-D) |
| `position_/state_mahalanobis_total_cov` | the same with C_obs + C_TLE,empirical. NaN unless an ensemble of ≥ 3 TLEs is available |

The distances are *squared* Mahalanobis distances, computed by Cholesky
decomposition. A distance is NaN when its covariance is not finite, not
positive definite, or has a condition number above 10¹⁴. The flag
`INVALID_STAGE5_COVARIANCE` and the precision-validation condition depend on
the 3×3 *position* block. If only the 6×6 matrix is too ill-conditioned (for
example because of r–v correlations close to 1 in a smoothed orbit), the
state D² alone is NaN. The run log notes this, and no flag is raised.

Stage-5 epochs are strongly correlated, so distances are reported per epoch
and summarised by medians. They are never summed into a χ² or p-value over
epochs.

## 7. Empirical reference dispersion (TLE ensemble)

The ensemble comes from `--tle-ensemble-files` and/or `--tle-ensemble-dir`
(`*.tle`, `*.txt`). It contains every distinct set of the same NORAD ID,
deduplicated by line content, with epoch within `--tle-ensemble-window-days`
(7) of the trajectory midpoint. The primary set is included when it falls in
the window, and the manifest says whether it did.

Each member is propagated and transformed exactly like the primary. At every
epoch, the sample mean and sample covariance (`ddof = 1`) of the 6-D states
give the empirical reference dispersion.

- With fewer than 3 members, `TLE_ENSEMBLE_SMALL` is raised and no
  reference-aware statistic is produced.
- With fewer than 7 members, the 6×6 sample covariance is rank-deficient. It
  is only ever used added to C_obs, which stays positive definite.
- Without a usable ensemble, `REFERENCE_UNCERTAINTY_UNKNOWN` is raised.
- `TLE_REFERENCE_DISPERSION_HIGH` marks rows where the ensemble's position
  dispersion exceeds the Stage-5 position σ.
- The primary TLE stays the reference. The ensemble mean never replaces it;
  the median |primary − mean| in position and velocity is reported
  (`summary.json` → `ensemble.primary_minus_ensemble_mean`).
- The full per-epoch 6×6 dispersion is stored in the residual table as
  `cov_tle_dispersion_<a>_<b>` (21 upper-triangle terms, GCRS). These are the
  position and velocity dispersion covariances, with their cross terms.

## 8. Summary statistics, trend and autocorrelation (`summary.json`)

For each RTN component, and for |Δr|, |Δv|, the angle, the range, the
along-track time and the distances, `summary.json` reports:

- mean, median and RMS;
- the MAD-based σ (1.4826·MAD);
- the 16th and 84th percentiles;
- the minimum and maximum.

**Trend.** `dr_j(t) = a_j + b_j (t − t0)` is fitted by weighted least
squares, with the Stage-5 RTN variance as weight. The slope
(`summary.statistics.<component>.trend.slope`, in m/s) is a mismatch *rate*,
not an orbit correction. Its formal σ assumes independent epochs and is for
reference only. `STRONG_RESIDUAL_TREND` is descriptive: the fitted change
over the arc exceeds 3× the median Stage-5 σ of that component.

**Autocorrelation.** Each mean-removed component gets a lag-1
autocorrelation and an effective sample size n(1 − ρ)/(1 + ρ).
`RESIDUAL_AUTOCORRELATED` uses the white-noise bound 3/√n. A smoothed
trajectory nearly always exceeds that bound, and for a pure offset or ramp ρ
is close to 1. The flag therefore marks structure; it is not evidence of a
problem, and it does not enter precision validation.

## 9. Precision validation

Every run states `precision_validation: true/false`. It is true only if all
of these hold:

1. the TLE was found and its NORAD ID matches;
2. the largest |age| over the compared epochs is ≤ `--tle-age-stale-days`;
3. the Stage-5 trajectory has at least one compared epoch;
4. the frame transformation succeeded;
5. no compared epoch has an unusable Stage-5 covariance.

`true` does not mean the TLE is exact. It means the comparison meets the
pipeline's minimum conditions for precision analysis. When it is false, every
row carries `NOT_PRECISION_VALIDATION`.

**Exit status:**

| Code | Meaning |
|---|---|
| 0 | precision validation true |
| 0 with `--allow-stale-reference` | a descriptive comparison, when TLE age is the *only* failed condition |
| 3 | completed but not valid for precision analysis |
| 2 | invalid input or configuration, including `TLE_NOT_FOUND` and `TLE_NORAD_MISMATCH` |
| 1 | unexpected failure, including `TLE_FRAME_FAILURE` (SGP4 or frame error) |

**Residual thresholds.** `HIGH_POSITION_RESIDUAL` and `HIGH_VELOCITY_RESIDUAL`
are set only when `--position-residual-warning-km` or
`--velocity-residual-warning-m-s` is given; there are no built-in kilometre
thresholds. `HIGH_CONDITIONAL_MAHALANOBIS` marks a position D² above
χ²₃(0.9973) = 14.16, conditional on the TLE.

## 10. Outputs

```text
results/stage6/<YYYYMMDD>/NORAD_<id>/
├── residuals.{parquet,csv}     one row per Stage-5 epoch (§6 columns plus observed and TLE states in GCRS and
│                               ECEF, RTN sigmas/covariances of observation and reference, flags)
├── tle_trajectory.parquet      TLE GCRS position and velocity, ECEF position, age, RTN basis vectors (for Stage 7)
├── summary.json                §8 statistics, TLE, ages, precision-validation conditions, ensemble, flag counts
├── stage6_manifest.json        inputs and hashes, TLE lines/file/SHA-256/source/epoch/age thresholds, SGP4
│                               version, frames, residual definition, policy, precision validation, git commit,
│                               command line, software versions, outputs
├── diagnostics/  rtn_residuals  position_residual_norm  velocity_residuals  angular_residual
│                 residual_vs_uncertainty  mahalanobis  tle_age  trajectory_overlay   (.png)
└── logs/stage6.log
```

A run can be reproduced from the Stage-5 manifest and trajectory, the exact
TLE (lines and SHA-256), the SGP4 and frame implementation, the git commit
and the command line. No network is needed.

**Light time.** Stages 4–5 do not correct light time. A measured GNSS
position therefore belongs to about t − 70 ms, which on its own gives an
along-track residual of about −0.27 km, or −0.07 s equivalent. The manifest
records this so that Stage 7 does not mistake it for a clock error.

## 11. Synthetic validation

`tests/stage6/` has 62 tests. The reference is the Stage-1 element set of
GPS BIIF-9 (PRN 26, NORAD 40534; epoch 2026-10-06T17:48Z). Observed
trajectories are built by the independent TEME route of §4, so a zero
residual is not circular.

| Specification test | Result |
|---|---|
| 34.1 zero residual | \|Δr\| < 0.1 m, \|Δv\| < 1 mm/s, through the full `tart-validate` path |
| 34.2 Cartesian offset | Δr equals the injected GCRS offset; `obs − tle` columns equal `dr_*` |
| 34.3 RTN offset (+1000, +2000, −500 m) | recovered to < 0.1 m |
| 34.4 velocity offset | recovered to < 1 mm/s. Time-shift signature: dr_T = \|v\|τ, dv_R = −GM/r²·τ |
| 34.5 covariance rotation | Monte Carlo of rotated samples; a covariance built in RTN is recovered exactly |
| 34.6 TLE age | fresh (0.1 d), warning (5 d), stale (20 d, exit 3), back-propagated (−1 d) |
| 34.7 single TLE | `position_mahalanobis_obs_cov` = 9 for a 3σ offset; `*_total_cov` NaN; `REFERENCE_UNCERTAINTY_UNKNOWN` |
| 34.8 ensemble | 5 nearby sets: dispersion mostly along-track; the total D² is below the conditional one; no unknown-reference flag |

The other tests cover:

- `--tle-file` selection from a history file;
- the Stage-1 provenance chain;
- an explicit download, with an injected fetcher and no network;
- the TLE-file failure cases: missing, wrong NORAD ID, modified file;
- mutually exclusive sources;
- the strict policy;
- user thresholds;
- plots;
- an unmodified Stage-5 run, checked by hashing every file before and after.

Synthetic nearby sets come from `sgp4init` plus `sgp4.exporter`. Epoch
shifts carry the mean anomaly, argument of perigee and RAAN forward with
SGP4's own secular rates, so an epoch-shifted set describes nearly the same
trajectory (0.5–1.8 km apart).

### Demonstration

Run `scripts/stage6_synthetic_demo.py`; outputs are in
[`examples/stage6`](examples/stage6), with the report in
[`demo_report.txt`](examples/stage6/demo_report.txt). The TLE epoch is
2026-10-06T17:48:11.977Z. Observations are 0.11 d after it, propagated
forward, so `precision_validation` is true and every run exits 0.

**Part A: known injection.** The observed trajectory is the TLE (by the
independent route) plus known offsets. The Stage-5-like covariance has σ
R/T/N = 300/800/400 m.

| Quantity | Value |
|---|---|
| injected RTN offset | (+1000.0, +2000.0, −500.0) m |
| recovered radial / along / cross | (+1000.000, +2000.001, −500.000) m |
| injected RTN velocity offset | (+0.100, −0.050, +0.020) m/s |
| recovered RTN velocity residual | (+0.1000, −0.0500, +0.0200) m/s |
| 3-D residual norm | 2291.3 m (expected 2291.3 m) |
| geocentric angular residual | 16.16″ |
| equivalent along-track time offset | 0.512 s |
| conditional position D² (TLE fixed) | 18.92, equal to (1000/300)² + (2000/800)² + (500/400)² |
| conditional state D² | 19.10 |
| single TLE | total-covariance D² NaN; `REFERENCE_UNCERTAINTY_UNKNOWN` |
| with a 5-set ensemble (primary included) | dispersion R/T/N = 25 / 2700 / 52 m; reference-aware position D² 13.31, state D² 14.73 |

**Part B: end to end through Stage 5.** The truth is the TLE trajectory
shifted by τ = 2 s, a pure along-track displacement of |v|τ = 7.818 km. It is
observed as Stage-4-like positions: 3 km along the line of sight, 300 m
across, 180 epochs over 30 min. These go through `tart-track` (two-body) and
then `tart-validate`, and both exit 0.

| Quantity | Value |
|---|---|
| recovered median R / T / N | (+0.080, +7.789, −0.017) km |
| median Stage-5 σ R / T / N | (0.272, 0.051, 0.051) km |
| equivalent along-track time offset | 1.992 s (injected 2 s) |
| velocity residual (median) | 1.178 m/s, of which dv_R = −1.177 m/s, matching −GM/r²·τ = −1.13 m/s |
| conditional position D² | median 6.7×10⁴: the time shift is a real discrepancy |

**Monte Carlo with τ = 0** (100 realisations × 5 epochs; no Stage-5
restarts). Here the TLE equals the truth, so the conditional position D²
should follow χ²₃.

| | Value | Expected (χ²₃) |
|---|---|---|
| mean / median | 2.82 / 2.33 | 3 / 2.37 |
| fraction ≤ 1 / 4 / 9 | 0.170 / 0.750 / 0.992 | 0.199 / 0.739 / 0.971 |

On white-noise synthetic data, Stage 5's uncertainties therefore carry
through Stage 6 calibrated. On the real data they do not (§12).

## 12. Real-data run (2026-06-25, GPS PRN 26, NORAD 40534)

Outputs are in [`examples/stage6_real`](examples/stage6_real).

**The only TLE available locally is not contemporaneous.** It is the Stage-1
set from the CelesTrak catalogue retrieved 2026-10-07: epoch 2026-10-06, about
**103 days after** the observation. CelesTrak serves only current elements,
and historical sets (for example from Space-Track) were not available
offline. Every run is therefore back-propagated 103 days:

- `TLE_STALE`, `TLE_BACKPROPAGATED` and `NOT_PRECISION_VALIDATION` are set;
- `precision_validation` is false.

These runs validate the software only. Their residuals must not be read as
TART astrometric accuracy, nor as TLE error.

| Run | Exit | Result |
|---|---|---|
| spec form, `--tle-file data/tle/gnss_2026-10-07T101156Z.tle`, two-station Stage-5 run (`diagnostic_min2`) | **3** | not precision validation (TLE age −103.35 d) |
| same, `--tle-source stage1 --allow-stale-reference` (`diagnostic_min2_stage1_allow_stale`) | 0 | descriptive. The Stage-1 provenance chain gives the same element set. |
| default Stage-5 run (`default`, no valid Stage-5 rows) | 3 | no epoch compared |

Descriptive residuals over the 36 compared epochs (09:18:35–09:24:22,
observed − TLE). These are not interpreted.

| Component | median | MAD σ | min … max | median Stage-5 σ | lag-1 ρ |
|---|---|---|---|---|---|
| radial | +2765 km | 285 km | +2400 … +3144 km | 129 km | 0.92 |
| along-track | +370 km | 33 km | +334 … +423 km | 8.6 km | 0.92 |
| cross-track | −627 km | 50 km | −693 … −563 km | 15.5 km | 0.92 |

- **3-D residual:** median 2860 km.
- **Geocentric angle:** 1.38–1.47°.
- **Equivalent along-track time:** 87–110 s.
- **Velocity residual:** median 2.32 km/s, about −2.2 km/s of it radial.
- **Conditional position D²:** median 3.1×10³.
- **Flags:** `HIGH_CONDITIONAL_MAHALANOBIS`, `STRONG_RESIDUAL_TREND` and
  `RESIDUAL_AUTOCORRELATED` everywhere, with `REFERENCE_UNCERTAINTY_UNKNOWN`
  (there is no ensemble).

What can be said without attribution:

- The residual is mostly radial, which matches the positive along-LOS offset
  Stage 4 found.
- It is hundreds of times the Stage-5 σ. Stage 5 already showed that its σ
  is optimistic on these data (`DYNAMICS_MISMATCH`).
- A 103-day back-propagated TLE is not a precision reference.

A scientific acceptance result needs an element set from close to
2026-06-25. Supply it with `--tle-file`, for example a Space-Track history
file for NORAD 40534; the nearest epoch is selected automatically. An
ensemble of the neighbouring sets in `--tle-ensemble-dir` would give the
reference-aware statistics.

## 13. Interpretation guardrails

A large residual does not by itself imply a bad TLE, a bad telescope, a
manoeuvre, a clock error, a calibration error or RFI. A small residual does
not prove absolute accuracy or correct uncertainty calibration. Stage 6
reports raw `observed − TLE` residuals and their uncertainty context;
systematic attribution is Stage 7.

## 14. Command name

Stage 6 adds `tart-validate` (`tart_gnss.stage6.cli:main`). The public
commands are:

- `tart-stations-update`
- `tart-stations-list`
- `tart-common-gnss`
- `tart-retrieve`
- `tart-localise`
- `tart-triangulate`
- `tart-track`
- `tart-validate`

There are no aliases. `tests/test_cli_names.py` checks that the entry points
match this set exactly, and that every `prog=` matches its entry point.
