# Stage 5: TLE-independent trajectory tracking (`tart-track`)

Stage 5 turns the epoch-by-epoch 3-D positions of Stage 4 into a continuous
satellite state: position and velocity with a full 6×6 covariance at every
epoch. It uses an extended Kalman filter (EKF) and a Rauch–Tung–Striebel
(RTS) smoother. The Stage-4 covariances enter the measurement update in full.

**No TLE enters Stage 5** — not as initial state, prior, process model,
regularisation target or tuning reference. Stage 4 copies a
`satellite_reference` block into its manifest for Stage 6, and Stage 5 never
reads it. A test makes every TLE loader and Skyfield's `EarthSatellite`
raise, then runs the full pipeline (§13).

Specification: [`tart_gnss_stage5_implementation.md`](../tart_gnss_stage5_implementation.md).

```text
Stage-4 positions + full 3x3 covariance (ITRS/ECEF)
        |  quality policy, covariance checks
        |  rotation to GCRS (Skyfield)
        v
TLE-independent initial state (first Stage-4 measurements)
        v
forward EKF (constant-velocity / two-body / J2), chi^2_3 innovation gate
        v
backward RTS smoother
        v
trajectory.parquet (GCRS state, ECEF position, 6x6 covariances) + innovations.parquet
```

## 1. Inputs

- `stage4_manifest.json`. It must have `stage == 4`, `schema_version == 1`
  and `coordinate_frame == "ITRS/ECEF"`, and its SHA-256 is recorded.
- `positions.{parquet,csv}`. Stage 5 looks next to the manifest first, then
  at the paths recorded in `outputs`, so a moved result directory still loads.

The required columns are those listed in §4 of the specification. The
measurement covariance is rebuilt in full from the six Stage-4 terms
(`cov_xx … cov_zz`). It is never reduced to its diagonal. §13 shows that
dropping the off-diagonal terms costs accuracy and calibration.

## 2. Frames

The state is propagated in GCRS, which is inertial. The measurements are in
ITRS/ECEF.

| Quantity | Transformation |
|---|---|
| rotation | `R(t) = skyfield.framelib.itrs.rotation_at(t)`, GCRS → ITRS: IAU 2000A/2006 precession–nutation, Greenwich apparent sidereal time (UT1), polar motion |
| position | `r_GCRS = R(t)ᵀ r_ITRS`; output `r_ITRS = R(t) r_GCRS` |
| covariance | `C_GCRS = R(t)ᵀ C_ITRS R(t)`. The transform is a pure rotation at fixed epoch, so this is exact |
| time | `load.timescale(builtin=True)`: Skyfield's bundled leap-second and UT1−UTC tables, so no network is needed. Unix seconds are split into whole days plus seconds of day, so leap seconds are handled. |

The builtin timescale has no polar-motion table, so polar motion is zero.
That costs at most about 0.3″, or about 40 m at GNSS range. UT1 comes from
Skyfield's bundled predictions. Both effects are negligible against
kilometre-scale Stage-4 uncertainties.

Only positions are transformed back to ECEF. An Earth-fixed velocity would
need the ω×r term, and it is not a Stage-5 output.

The frame code is checked in three ways:

- `M r_ITRS` matches Skyfield's own GCRS position of a ground point,
  `wgs84.latlon(...).at(t)`, to within 1 mm.
- The matrices are rotations, and eigenvalues survive the round trip.
- Rotated noise samples have covariance `M C Mᵀ`.

## 3. Time handling

- The real Stage-4 timestamps are used, with the actual Δt at every step.
  Uniform cadence is never assumed.
- Rows out of time order are sorted, and the log and manifest say so
  (`input_time_ordered`).
- For an exactly repeated timestamp, the first row is used and later ones are
  rejected as `DUPLICATE_TIMESTAMP`.
- `LARGE_TIME_GAP` marks a row if more than `--max-gap-s` (default 120 s)
  has passed since the last assimilated measurement.
- Gaps need no special treatment. The process noise grows as Δt²–Δt⁴ and the
  dynamics propagate the covariance, so prediction uncertainty widens on its
  own (§13).

## 4. Dynamical models

The state is `x = [x, y, z, vx, vy, vz]`, in metres and m/s, in GCRS.

| `--dynamics` | Model | Transition |
|---|---|---|
| `constant-velocity` | r(t+Δt) = r + vΔt, v constant | exact, `F = [[I, ΔtI], [0, I]]` (diagnostic baseline) |
| `two-body` (default) | r̈ = −GM r/\|r\|³ | `solve_ivp` (DOP853, rtol 1e-10, atol 1e-6) on the state plus the variational equations Φ̇ = AΦ, `A = [[0, I], [G, 0]]`, `G = GM/r³ (3 r̂r̂ᵀ − I)` |
| `j2` (optional) | two-body plus J2 about the true rotation pole at the first epoch (ITRS z axis in GCRS) | as above; the J2 part of `G` comes from central differences of the analytic J2 acceleration |

The constants come from IERS Conventions (2010), Table 1.1:

- GM = 3.986004418×10¹⁴ m³ s⁻², which equals the WGS84 value;
- J2 = 1.0826359×10⁻³;
- a_E = 6 378 136.6 m, which is identical to `skyfield.constants.ERAD`.

The tests cross-check GM against `sgp4.earth_gravity.wgs84.mu`, used as a
constant only.

Because Φ is integrated with the state, it is the exact Jacobian of the
integrated flow. The integrator is validated three ways:

- against an *analytic* Kepler solution (Lagrange f and g with the universal
  Kepler equation): agreement to 1 mm over 1 h;
- Φ against central differences of that analytic solution: relative error
  below 10⁻⁷;
- conservation of energy and angular momentum over 2 h, to 10⁻¹¹.

The synthetic truths are generated with the same analytic solution, so the
filter's integrator never checks itself.

## 5. Process noise

`--process-accel-sigma` (σₐ, m/s²) is white acceleration noise held constant
over each interval. This is the form given in the specification:

```text
Q(Δt) = σₐ² [[Δt⁴/4 I, Δt³/2 I], [Δt³/2 I, Δt² I]]
```

The orbital models use the same Q, to first order: over seconds to minutes
the gravity gradient changes Q only at second order in Δt.

The defaults are set by the size of the accelerations each model leaves out.
They are never tuned against TLE agreement.

| Model | Default σₐ | Reason |
|---|---|---|
| constant-velocity | 1 m/s² | gravity itself is unmodelled (GM/r² ≈ 0.57 m/s² at GNSS radius) |
| two-body | 1e-4 m/s² | J2 ≈ 5e-5 m/s² at GNSS radius, lunisolar ≈ 5e-6, radiation pressure ≈ 1e-7 |
| j2 | 1e-5 m/s² | lunisolar terms dominate what is left |

The value used and its source (option or model default) go in the manifest.

## 6. TLE-independent initialisation

Only Stage-4 measurements are used. `t0` is the first usable measurement.
The measurements within `--initial-velocity-window` (default 30 s) of `t0`
are fitted by generalised least squares (GLS) with their full covariances:

```text
y_k = r0 + v0 (t_k − t0) + ½ a(r0) (t_k − t0)²
```

Here `a` is the model acceleration, zero for constant velocity. It is a
known correction, not a fitted parameter. The fit gives (r0, v0) and their
6×6 covariance N⁻¹.

The method falls back in steps:

1. ≥ 3 points: `linear_fit`.
2. 2 points (from the window, else the first two measurements):
   `finite_difference`, which is the GLS fit through two points.
3. 1 point: `zero_velocity`, with 10 km/s σ per axis.

`INITIAL_VELOCITY_WEAK` is raised for anything other than a `linear_fit`, or
when any velocity σ exceeds 1 km/s.

**Avoiding double counting.** The filter assimilates the window measurements
again. If the fit covariance N⁻¹ were also the prior, those data would count
twice and the early states would be over-confident. The prior is therefore
`P0 = 10⁴ N⁻¹`. The fit sets the linearisation point and the *shape* of the
initial uncertainty, including the r–v correlation, but adds at most 10⁻⁴ of
the window's information.

The early-epoch filtered NEES checks this, as it is where double counting
would appear:

- 300 runs: NEES at the first epochs is 5.7–6.2, with χ²₆ mean 6.
- At epoch 0 it is 2.95, because the velocity there is still effectively
  unconstrained, leaving 3 effective degrees of freedom.

`tests/stage5/test_covariance_consistency.py` asserts it.

## 7. Extended Kalman filter

The measurement model is `y_k = H x_k + v_k` with `H = [I₃ 0]`. `R_k` is the
full Stage-4 covariance rotated to GCRS.

```text
prediction   x_k|k-1 = f(x_k-1|k-1, Δt_k),  P_k|k-1 = F_k P F_kᵀ + Q(Δt_k)
innovation   ν_k = y_k − H x_k|k-1,          S_k = H P_k|k-1 Hᵀ + R_k
NIS          D²_k = ν_kᵀ S_k⁻¹ ν_k  ~ χ²₃
gain         K_k = P_k|k-1 Hᵀ S_k⁻¹         (Cholesky solves, no explicit inverse)
update       x_k|k = x_k|k-1 + K_k ν_k
             P_k|k = (I − KH) P (I − KH)ᵀ + K R Kᵀ    (Joseph form)
```

After every prediction and update, P is symmetrised and checked: it must be
finite, with smallest eigenvalue ≥ −10⁻⁹ × the largest. A failure raises
`FilterDivergenceError` (exit 3) instead of letting the covariance collapse
silently.

**Gating** (`--innovation-gate-probability`, default 0.9973, giving a
threshold of χ²₃ = 14.16):

- A measurement above the gate is an `INNOVATION_OUTLIER`. It stays in both
  output tables.
- With `--robust-update none` (the default) it is not assimilated.
- With `inflate` it is assimilated with `R` scaled by `D²/gate`, which brings
  its NIS down to the gate.

**Divergence recovery.** With near-deterministic two-body dynamics, a single
early gate decision can lock the filter onto a wrong state. The predicted
error and its σ then grow together, and every later measurement is gated
too.

The 200-realisation Monte Carlo of the demo found this once. A
0.1%-probability fourth measurement (NIS 15.7) was gated, and the run ended
2 400 km off.

The fix: after 3 consecutive gated measurements, the filter restarts at the
first of them, re-initialised from the data (§6).

- Those rows are flagged `FILTER_DIVERGENCE` (advisory, because the restarted
  filter assimilates them).
- The restart is a segment boundary, and the smoother does not cross it.
- The manifest records `n_resets` and the reset times.
- The consistency statistics use each measurement's *first*, sequential NIS,
  so a restart cannot hide the evidence that triggered it (§9).

## 8. RTS smoother

```text
C_k   = P_k|k F_k+1ᵀ P_k+1|k⁻¹     (linear solve)
x_k|N = x_k|k + C_k (x_k+1|N − x_k+1|k)
P_k|N = P_k|k + C_k (P_k+1|N − P_k+1|k) C_kᵀ
```

`x_k+1|k` is the stored *nonlinear* prediction `f(x_k|k)`, and `F_k+1` is
the flow Jacobian.

The smoothed state is the primary offline product. `--no-smoother` reports
filtered states only. If the smoother fails, every row gets
`SMOOTHER_FAILED`, which is blocking.

In the linear case (constant velocity, no process noise), the smoothed state
at t0 equals the batch GLS fit of all the data. The tests check this, and
check that the smoothed covariance never exceeds the filtered one.

## 9. Innovation statistics (TLE-independent consistency)

The statistics cover the measurements offered by the quality policy, gated
ones included. They use each measurement's first NIS. The first 2 innovations
are skipped as burn-in, because they fix the six state parameters from an
almost uninformative prior and have NIS ≈ 0.

| Statistic | Expected for a consistent model |
|---|---|
| mean / median NIS | 3 / 2.37 (χ²₃) |
| fraction above the gate | 1 − gate probability |
| lag-j autocorrelation of whitened innovations `L⁻¹ν` (`S = LLᵀ`) | 0 ± 1/√(3 n_pairs) |

Run-level flags, which are applied to every row:

- **`HIGH_NIS`** is raised when either test fails:
  - the mean NIS *inside the gate* exceeds the 0.9973 quantile of χ²₃ₙ/n;
  - more innovations exceed the gate than the 0.9973 quantile of
    Binomial(n, 0.0027) allows.

  One gated outlier therefore does not raise it. The reported `mean_nis`
  still includes such outliers.
- **`DYNAMICS_MISMATCH`** is raised when the lag-1 autocorrelation exceeds
  3/√(3 n_pairs). Correlated innovations mean unmodelled dynamics *or
  time-correlated measurement errors*.

## 10. Stage-4 quality policy and Stage-5 flags

| Stage-4 condition | `strict` (default) | `permissive` |
|---|---|---|
| `TRIANGULATION_FAILED`, `INSUFFICIENT_STATIONS`, `NEGATIVE_RANGE`, `COVARIANCE_SINGULAR`, `fit_success` false, non-finite position, invalid covariance, duplicate time | rejected | rejected |
| `POOR_GEOMETRY` (Stage-4 blocking), `STAGE3_QUALITY_WARNING` (Stage-4 epoch built from Stage-3 rows Stage 3 marked unusable) | rejected | used, `STAGE4_QUALITY_WARNING` |
| `TIME_MISMATCH`, `HIGH_LINE_RESIDUAL`, `HIGH_CHI2`, `INPUT_DIRECTION_WARNING`, `CALIBRATION_WARNING`, `BOOTSTRAP_UNSTABLE` | used, `STAGE4_QUALITY_WARNING` | used, `STAGE4_QUALITY_WARNING` |

`strict` therefore uses exactly the Stage-4 `valid` epochs, minus those built
from unusable Stage-3 rows. On 2026-06-25, those come from a different source
([Stage-4 doc §12](stage4_triangulation.md)).

**Covariance integrity.** Each Stage-4 covariance must pass four checks:

- it is finite;
- it is symmetric to a relative tolerance of 10⁻⁶ (smaller asymmetry is
  removed by `(R + Rᵀ)/2`);
- it is positive definite;
- its condition number is ≤ 10¹².

A failing row is rejected and flagged `INVALID_STAGE4_COVARIANCE`. A
covariance is never repaired.

| Stage-5 flag | Meaning | Blocking |
|---|---|---|
| `INSUFFICIENT_MEASUREMENTS` | fewer than 3 usable measurements (no filter run) | yes |
| `NOT_INITIALISED` | before the first usable measurement (no state) | yes |
| `SMOOTHER_FAILED` | RTS smoother failed | yes |
| `COVARIANCE_INVALID` | reserved: a state-covariance failure currently aborts the run (`FilterDivergenceError`, exit 3) | yes |
| `FILTER_DIVERGENCE` | one of ≥ 3 consecutive gated measurements that triggered a restart | no |
| `INVALID_STAGE4_COVARIANCE` | the Stage-4 covariance failed the checks | no |
| `MEASUREMENT_REJECTED` / `NO_MEASUREMENT` | row not used by the policy / no finite Stage-4 position | no |
| `INNOVATION_OUTLIER` | NIS above the gate | no |
| `LARGE_TIME_GAP` | more than `--max-gap-s` since the last assimilated measurement | no |
| `EXTRAPOLATED` | after the last assimilated measurement (prediction only) | no |
| `INITIAL_VELOCITY_WEAK` | within the data behind a weak initial velocity | no |
| `HIGH_NIS`, `DYNAMICS_MISMATCH` | run-level consistency failures (§9) | no |
| `STAGE4_QUALITY_WARNING` | Stage-4 row used with advisory flags | no |
| `GOOD` | none of the above | — |

**Exit status:**

| Code | Meaning |
|---|---|
| 0 | ≥ 3 measurements assimilated and ≥ 1 valid row |
| 3 | completed but not enough valid trajectory estimates (including filter divergence) |
| 2 | invalid input or configuration |
| 1 | unexpected failure |

## 11. Outputs

```text
results/stage5/<YYYYMMDD>/NORAD_<id>/
├── trajectory.{parquet,csv}     one row per Stage-4 epoch
├── innovations.{parquet,csv}    one row per epoch with a finite Stage-4 position
├── stage5_manifest.json
├── diagnostics/  position_track  velocity_track  uncertainty_track  innovations
│                 nis  stage4_vs_smoothed  covariance_eigenvalues   (.png)
└── logs/stage5.log
```

**`trajectory`** holds:

- `time_utc` and `time_mjd`;
- the filtered and smoothed GCRS position and velocity;
- the smoothed ECEF position (rotated back);
- `speed_smoothed_m_s`;
- the filtered and smoothed position σ, and the smoothed velocity σ (all
  GCRS);
- `acceleration_m_s2`, the model |a|, diagnostic only: GM/r² for two-body, 0
  for constant velocity;
- the Stage-4 position in ECEF and GCRS, and its GCRS σ;
- `measurement_used`, `measurement_accepted`, `rejection_reason`,
  `time_since_update_s`, `nis`, `valid`, `quality_flags` and the Stage-4
  flags;
- the smoothed 3×3 ECEF position covariance;
- **the full filtered and smoothed 6×6 covariances**, as 21 upper-triangle
  columns each: `cov_{filtered,smoothed}_<a>_<b>`, with a, b in
  x y z vx vy vz. `tart_gnss.stage5.output.state_covariance(row, kind)`
  rebuilds the matrix.

**`innovations`** has one row per epoch with a finite Stage-4 position,
from the first usable epoch on (earlier epochs have no prediction). It holds,
per row:

- ν and S (six terms);
- the whitened innovation;
- `nis` and `nis_first_pass` (they differ only for rows re-processed after a
  restart);
- the gate, `accepted`, `attempted` and `inflation`;
- the Stage-4 and Stage-5 flags.

**`stage5_manifest.json`** holds:

- the exact inputs and their SHA-256;
- the frames and the frame implementation;
- `tle_used: false`;
- the state model: constants and their source, integrator and tolerances,
  STM method, J2 pole;
- the filter: initialisation including fit σ, prior inflation, divergence
  recovery and resets;
- the smoother, process noise, innovation gate, innovation statistics and
  dynamics diagnostics;
- the quality policy with its rejection counts;
- a summary with flag counts;
- Stage-4 provenance;
- software versions with the git commit and dirty state;
- the command line and output paths.

## 12. Stage-6 interface

At every epoch Stage 6 can read:

- `time_utc` / `time_mjd`;
- the smoothed GCRS position and velocity;
- the smoothed ECEF position, with its 3×3 covariance;
- the full 6×6 smoothed state covariance (GCRS);
- `valid`, `quality_flags`, `nis` and the innovations table.

Stage 6 will propagate the archived TLE to these epochs on its own. Stage 5
computes no TLE residual.

## 13. Synthetic validation

**Truth.** A circular GNSS-like orbit: radius 26 560 km, inclination 55°,
above the example stations at 2026-06-25T09:14Z. It is propagated with the
analytic Kepler solution.

**Stage-4-like measurements.** The truth plus Gaussian noise from an
anisotropic covariance: σ_LOS along the line of sight from the South-Africa
station, σ_T across it. This is the shape of real Stage-4 covariances.

**Tests** (`tests/stage5/`, 90 in total). Each line below is one of the
specification's §32 tests, or a related check:

- **32.1 Exact two-body.** 60 epochs at 5 km / 1 km. Smoothed position RMSE
  < 1.5 km, below 3× the predicted σ. Velocity RMSE < 3× the predicted σ.
- **32.2 Anisotropy.** σ_LOS 50 km, σ_T 1 km, 20 realisations.
  - The full covariance gives 0.77× the position error of the diagonal-only
    approximation.
  - The full covariance is calibrated (position NEES 3.3 against 3). The
    diagonal approximation is not (4.7).
  - The smoothed covariance's long axis lies along the line of sight
    (|cos| > 0.95).
- **32.3 Irregular cadence.** Random epochs: the true Δt is used and NIS is
  consistent. Labelling jittered data with a nominal grid raises mean NIS
  more than 10×.
- **32.4 Missing epochs.** A 200 s gap:
  - `LARGE_TIME_GAP` is raised, with `time_since_update_s` = 210 s;
  - the predicted covariance across the gap is more than 4× the last
    filtered one (trace);
  - the filtered σ grows monotonically through the gap;
  - the smoothed σ in the gap exceeds that of the same epoch with full data.
- **32.5 Outlier.** One +80 km measurement is the only one gated. It is kept
  in both tables, flagged and still `valid`. `--robust-update inflate`
  assimilates it with inflated R.
- **32.6 Process-noise mismatch.** The truth has an unmodelled constant
  5×10⁻³ m/s² acceleration.
  - With σₐ = 10⁻⁷, gate exceedances are far above chance, the filter has to
    restart, and `HIGH_NIS` + `DYNAMICS_MISMATCH` are raised.
  - With σₐ = 5×10⁻³ there are no exceedances or restarts, and both NIS and
    autocorrelation are lower.
  - A *constant* acceleration is not white noise, so some autocorrelation
    remains even with matched σₐ. That is correct behaviour.
- **32.7 No TLE.** Every TLE loader and `EarthSatellite.__init__` raises, and
  the Stage-4 manifest carries a TLE block. The pipeline still runs, exits 0,
  and records `tle_used: false`.
- **Two-body against constant velocity**, both with σₐ = 10⁻⁴, 1 km / 300 m
  noise at 20 s cadence. Two-body has:
  - under half the innovation RMS;
  - lower median NIS;
  - no lower an accepted fraction;
  - no autocorrelation, where constant velocity shows ≈ 0.56;
  - no restarts, where constant velocity needs 4.
- **Lock-out regression.** A 300 km jump in the measured track triggers
  exactly one restart, at the jump. Both sides of it are tracked within
  30 km.
- **Monte Carlo** (`test_covariance_consistency.py`, 60 runs): mean smoothed
  NEES and early filtered NEES are within 3.5 sd of 6, and the 1/2/3σ
  coverage of position and velocity matches χ²₃.

### Demonstration

Run `scripts/stage5_synthetic_demo.py`; outputs are in
[`examples/stage5`](examples/stage5), with the report in
[`demo_report.txt`](examples/stage5/demo_report.txt).

**Part A** uses 156 epochs over 30 min:

- a jittered 10 s cadence;
- a 240 s gap;
- one +150 km outlier;
- σ_LOS = 30 km, σ_T = 3 km, so the median single-epoch 3-D error is
  20.4 km.

| Quantity | Value |
|---|---|
| true initial position (GCRS, km) | (345.728, 15930.345, −21249.428) |
| true initial velocity (GCRS, m/s) | (−3683.21, 988.64, 681.24) |
| smoothed initial position (km) | (345.884, 15927.914, −21244.362), off by 5.6 km |
| smoothed initial velocity (m/s) | (−3682.76, 988.28, 680.21), off by 1.18 m/s |
| position / velocity RMSE, all epochs | 5.18 km / 1.12 m/s |
| predicted smoothed 3-D σ (median) | 3.06 km (max 4.1 km, in the gap); 3.6 m/s |
| mean / median NIS | 15.45 / 2.44. The mean includes the gated 150 km outlier, and `HIGH_NIS` is correctly not raised |
| innovation rejection | 2 of 156: the injected outlier (NIS 1891) and one chance exceedance at epoch 4 (NIS 17.8) |
| lag-1 autocorrelation | −0.011 (limit 0.140) |

**Monte Carlo** (200 realisations; 4 epochs each; no outlier).

One realisation needed a divergence restart. By selection, its pre-restart
segment is a bad draw, so it is excluded below; with it, mean NEES is 6.63.

| | Empirical | Predicted / expected |
|---|---|---|
| 3-D position error (RMS) | 3.53 km | 3.45 km |
| 3-D velocity error (RMS) | 3.62 m/s | 3.61 m/s |
| mean / median NEES | 6.61 / 5.81 | 6 / 5.35 (χ²₆) |
| position coverage 1σ / 2σ / 3σ | 0.182 / 0.702 / 0.962 | 0.199 / 0.739 / 0.971 (χ²₃) |
| velocity coverage 1σ / 2σ / 3σ | 0.178 / 0.717 / 0.950 | 0.199 / 0.739 / 0.971 |

The epochs within one run are strongly correlated, so batch means scatter
widely. Two further 100-run batches of the same setup (every 10th epoch,
restarted runs excluded) gave mean NEES 6.24 and 5.73. The 6.61 above is
therefore within about 2.5σ of 6, given the batch scatter. A scratch check
with variants of the setup gave 5.25–6.29 (σₐ = 0, isotropic 10 km noise, no
gap); that script is not committed.
The predicted σ matches the empirical error to within about 2% for position and
0.3% for velocity.

**Dynamics comparison** (same data; TLE-independent diagnostics only):

| Model | σₐ (m/s²) | innovation RMS | mean / median NIS | lag-1 acf | accepted | position RMSE |
|---|---|---|---|---|---|---|
| two-body | 1e-4 | 35.40 km | 15.45 / 2.44 | −0.011 | 0.987 | 5.18 km |
| constant velocity | 1 (default) | 37.55 km | 13.92 / 2.54 | 0.010 | 0.987 | 8.80 km |
| constant velocity | 1e-4 | 43.25 km | 16.61 / 3.14 | 0.042 | 0.981 | 15.40 km |

At 30 km noise, constant velocity with its large default σₐ is *consistent*,
but less precise: 1.7× the position error. Constant velocity with the
two-body σₐ is worse on every diagnostic and has 3× the position error.

**Part B (end to end).** The same orbit is observed for 15 min (90 epochs,
elevation ≥ 43°) as directions from the three example stations, with 0.02°
noise. It is then run through `tart-triangulate` and `tart-track`, and both
exit 0.

| Quantity | Value |
|---|---|
| Stage-4 single-epoch 3-D error (median) | 39.0 km |
| Stage-5 smoothed position / velocity RMSE | 6.92 km / 10.50 m/s |
| predicted 3-D σ (median) | 8.45 km |
| mean / median NIS | 2.81 / 1.92 |
| lag-1 autocorrelation | 0.033 |
| run flags | none |

## 14. Real-data run (2026-06-25, GPS PRN 26, NORAD 40534)

Outputs are in [`examples/stage5_real`](examples/stage5_real). The inputs are
the Stage-4 example runs in [`examples/stage4_real`](examples/stage4_real):
`results/stage4/…` from the spec command is not kept in the repository, and
those runs are its committed copy.

**Specification command, on the default Stage-4 run.** This exits **3**, with
`INSUFFICIENT_MEASUREMENTS`. That run has a single valid three-station epoch
(09:23:54.9); the other 68 epochs are `INSUFFICIENT_STATIONS`. One position
cannot give a trajectory, and Stage 5 says so rather than inventing one.

**Two-station diagnostic Stage-4 run** (`--min-stations 2`). This exits 0.
It has 35 Stage-4 epochs from 09:18:35 to 09:24:15 (340 s, 10 s cadence),
with 34 assimilated.

| Quantity | Two-body | Constant velocity |
|---|---|---|
| mean / median NIS | 2.78 / 1.97 | 2.71 / 1.96 |
| lag-1 / 2 / 3 autocorrelation (limit 0.31) | 0.69 / 0.60 / 0.51 | 0.69 / 0.59 / 0.51 |
| innovation RMS | 248.6 km | 248.6 km |
| gated | 1 of 35 | 1 of 35 |
| median smoothed σ (GCRS x, y, z) | 8, 71, 106 km | 8, 71, 106 km |
| Stage-4 minus smoothed (RMS) | 269 km | 268 km |
| run flags | `DYNAMICS_MISMATCH` | `DYNAMICS_MISMATCH` |

What this shows:

- **The epoch-to-epoch scatter matches the Stage-4 covariance, but the errors
  are not white.** NIS is close to χ²₃. NIS cannot see an error common to
  all epochs, such as the ~3 000 km along-LOS offset found in Stage 4. However, the innovations are strongly
  autocorrelated out to lag 3 (`DYNAMICS_MISMATCH`). The Stage-4 positions
  wander about the smoothed track in a slow, correlated pattern
  (`stage4_vs_smoothed.png`, `position_track.png`). This matches the
  slowly varying 0.1–0.2° Stage-3 direction systematics. The filter assumes
  independent measurement errors, so its smoothed covariance is too
  optimistic here. A 340 s average of correlated errors does not shrink as
  1/√N.
- **The gated measurement is the run's only three-station epoch**
  (09:23:54.9, NIS 36). Its Stage-4 solution already carried
  `HIGH_LINE_RESIDUAL` from MG-Tana's line. Stage 5 independently finds it
  inconsistent with the two-station track.
- **Two-body against constant velocity is inconclusive.** At about 250 km
  per-epoch noise over 340 s, the curvature from gravity (at most ½·0.45 m/s²·340² ≈ 26 km) is far below
  the noise. The two models give the same diagnostics, so the default stays
  two-body on physical grounds, not because the data prefer it.
- **A TLE-independent plausibility check.** The smoothed radius is
  29 236–29 974 km, and the speed is 4.99–5.05 km/s ± 0.40 km/s (1σ along
  the velocity). A bound circular orbit at that radius moves at 3.67 km/s,
  and escape speed is 5.19 km/s. A speed 3.4σ above circular is further
  evidence of correlated Stage-4 errors that the formal covariance does not
  capture. This check uses only the estimated state and GM, not a TLE.
- The initial velocity comes from 4 points over 30 s, with σ = (2.1, 17.4,
  25.0) km/s. It is `INITIAL_VELOCITY_WEAK`, as expected at this noise
  level. The smoother removes the dependence on it.

No TLE residuals are computed here; that is Stage 6.

## 15. Limitations

- **The filter assumes white measurement errors.** Correlated Stage-4 errors
  make the smoothed covariance optimistic. Stage 5 detects this
  (`DYNAMICS_MISMATCH`) but does not model it. That would need a
  measurement-bias state (Stage 7), which is out of scope here.
- Process noise is first order for the orbital models. J2 uses a fixed pole
  for the whole arc, which is a negligible approximation over minutes.
- Positions only are rotated to ECEF; there is no ECEF velocity.
- Light time and aberration are not modelled, in Stage 4 or here. For GNSS
  they amount to ~70 ms, or about 0.3 km along track.
- A divergence restart splits the smoothed track into independent segments.
  Each is consistent on its own, but the segment that ends at a restart is,
  by selection, usually a poor estimate.

## 16. Command name

Stage 5 adds `tart-track` (`tart_gnss.stage5.cli:main`). The public commands
are:

- `tart-stations-update`
- `tart-stations-list`
- `tart-common-gnss`
- `tart-retrieve`
- `tart-localise`
- `tart-triangulate`
- `tart-track`

There are no aliases. `tests/test_cli_names.py` checks that the entry points
match this set exactly, and that every `prog=` matches its entry point.
