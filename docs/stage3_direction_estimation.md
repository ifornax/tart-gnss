# Stage 3: visibility-domain direction estimation

`tart-localise` measures the apparent direction of a GNSS satellite at each
TART station, independently, from the complex visibilities in the Stage-2
Measurement Sets. The primary estimator fits the visibilities directly. A
dirty-image peak finder is an independent cross-check. The outputs are
time-resolved (Az, El) tracks with covariances, quality metrics and flags,
which are the input to Stage 4 (multi-station localisation).

Stations are never combined into one interferometer, and their timestamps are
never forced onto a common grid. Time alignment belongs to Stage 4.

Specification: [`tart_gnss_stage3_implementation.md`](../tart_gnss_stage3_implementation.md).
Code: [`src/tart_gnss/stage3/`](../src/tart_gnss/stage3).

---

## 1. Scientific objective

For each station and each time bin, estimate the unit vector **s** toward the
dominant point-like source near the target satellite, with an uncertainty, from
the measured visibilities alone. The TLE prediction is used only to:

- choose the target satellite;
- centre a broad search region (default radius 10°);
- after the fit, flag detections of obviously unrelated sources
  (`TLE_UNASSOCIATED`, more than 5° away);
- plot a reference track.

It never enters the objective function. It is never a prior, never replaces a
failed fit, and is never written out as an observed position. The columns
`az_tle_deg`, `el_tle_deg`, `range_tle_km` and `tle_separation_deg` are
labelled reference predictions.

## 2. Measurement equation

For antennas *p* = `ANTENNA1` and *q* = `ANTENNA2` at ENU positions
**r**<sub>p</sub>, **r**<sub>q</sub>, define the baseline
**b**<sub>pq</sub> = **r**<sub>q</sub> − **r**<sub>p</sub> (metres). A single far-field
point source in direction **s** gives

```
V_pq(ν) = C · exp(−i φ_pq) + ε_pq,      φ_pq = 2π (ν / c) b_pq · s
```

Here *C* = *a* + *ib* is a complex nuisance amplitude, absorbing the source
power and any residual common phase. ν is each sample's own channel
frequency, so nothing assumes a single channel.

### 2.1 UVW and sign convention (validated)

`tart2ms` writes one FIELD per timestamp, with its phase centre at the local
zenith. It stores

```
UVW = r_p − r_q = −b_pq        (w ≈ 0 for the zenith phase centre)
```

so, in terms of the stored UVW, the same equation reads
`V = C exp(+2πi (u l + v m + w n)/λ)`. This is the opposite sign to the
textbook `exp(−2πi …)` form if UVW were plugged in naively. The code
therefore never uses the UVW column in the model. ENU baselines come from
the ANTENNA table positions, rotated to ENU at the station's geodetic
coordinates. The loader checks for every MS that the stored UVW equals −**b**
within 0.02 m + 2 % of the longest baseline, and that |w| is within the same
tolerance. If either check fails, the MS is rejected with
`InvalidVisibilityDataError`.

**Why a synthetic test cannot settle this.** A synthetic round trip uses the
same forward model as the estimator, so it recovers its own convention
whatever that is. A sign error maps (Az, El) to (Az + 180°, El), which looks
plausible. The convention was established and is regression-tested
(`tests/stage3/test_real_data.py`, run when the data are present) on the real
2026-06-25 Measurement Sets:

| Check | Brightest dirty-image peak (whole sky) | TLE prediction | Mirror (wrong sign) |
|---|---|---|---|
| ZA-Rhodes, 09:13:08.7 UTC | Az 175°, El 70° | GPS PRN 26 (target): Az 176.5°, El 71.0° | Az 355° (nothing predicted) |
| MU-UDM, 09:18:59 UTC | Az 110°, El 67° | GPS PRN 28: Az 106.2°, El 66.6° | Az 290° (nothing predicted) |
| MG-Tana, 09:19:17 UTC | Az 108°, El 55° | GPS PRN 28: Az 106.0°, El 54.4° | Az 288° |

On every station, the stored UVW matches −**b** to ≤ 0.016 m, with
|w| < 3 × 10⁻⁵ m.

## 3. Coordinate conventions

Local frame: x = East, y = North, z = Up (ENU). A direction is
**s** = (*l*, *m*, *n*), with *n* = √(1 − *l*² − *m*²) ≥ 0. *l* and *m* are
the East and North direction cosines (a zenith-centred tangent plane).

```
l = sin Az cos El,   m = cos Az cos El,   n = sin El
Az = atan2(l, m)  (0° = North, 90° = East),   El = arcsin n  (0° = horizon, 90° = zenith)
```

ECEF (ITRS) → ENU at geodetic latitude φ and longitude λ:

```
[E]   [ −sin λ          cos λ           0    ] [ΔX]
[N] = [ −sin φ cos λ   −sin φ sin λ    cos φ ] [ΔY]
[U]   [  cos φ cos λ    cos φ sin λ    sin φ ] [ΔZ]
```

Δ is taken relative to the mean antenna position. Unit tests cover the
cardinal directions, round trips, the up axis, and a blind-search check that
a source in the East is not recovered in the West (and likewise for North and
South).

## 4. Estimator

### 4.1 Complex-amplitude elimination

For a fixed direction, the weighted least-squares amplitude and the
profiled objective are

```
Ĉ(s)   = Σ w_k V_k e^{+iφ_k} / Σ w_k
χ²(s)  = Σ w_k |V_k − Ĉ e^{−iφ_k}|² = Σ w_k |V_k|² − |Σ w_k V_k e^{+iφ_k}|² / Σ w_k
```

which leaves a two-dimensional problem in (*l*, *m*). The complex visibility
is fitted directly. Wrapped phases are never subtracted, so there are no
phase-wrap discontinuities.

### 4.2 Coarse search

χ² is evaluated on a grid uniformly spaced on the sphere: exact great-circle
offsets on a square tangent-plane lattice about the cap centre. The grid step
is λ<sub>min</sub>/B<sub>max</sub>/5 (about 0.4–0.6° for TART), covering a
cap of radius `--search-radius-deg` (default 10°) about the TLE direction at
that bin's time. `--blind` searches the whole visible hemisphere instead.

### 4.3 Local refinement

All four parameters (*l*, *m*, Re *C*, Im *C*) are refined with
`scipy.optimize.least_squares` (`trf`, bounds −1 ≤ *l*, *m* ≤ 1), starting
from the best grid point and its Ĉ. The fit uses the stacked real residual
vector [Re *r*, Im *r*], *r* = √w (V − C e<sup>−iφ</sup>), with an analytic
Jacobian. *l*² + *m*² ≤ 1 is enforced by radial clipping inside the residual
function, and a solution on the unit circle is flagged.

### 4.4 Robust loss

`--robust-loss linear | huber | soft_l1` (default `huber`). For the robust
losses the amplitude is refined jointly: Ĉ in §4.1 minimises only the linear
loss. The loss scale is set from the data:
`f_scale = 1.345 × σ_robust`, where σ<sub>robust</sub> = 1.4826 MAD of the
residual components of the coarse solution. Without this, TART visibilities
(|V| ≪ 1) would make Huber identical to least squares. The loss and
`f_scale` are recorded. In the synthetic test with 10 % corrupted samples,
Huber and soft-L1 errors are less than half the least-squares errors.

### 4.5 Search boundaries

`BOUNDARY_SOLUTION` is set if any of the following holds:

- the refined solution is within half a grid step of the cap edge;
- the coarse optimum lies on the edge;
- χ² in a guard band two resolution elements wide just outside the cap is
  lower than anywhere inside.

In each case the true optimum probably lies outside the cap, and the in-cap
solution may be a sidelobe. Such solutions are reported but marked not
usable. They are never replaced by the cap centre. Sources further outside
than the guard band cannot be detected by a bounded search.

## 5. Multi-frequency handling

Each (row, channel) sample is a separate term with its own ν in φ. No
channel averaging is done. TART MSs have one 2.5 MHz channel, and
bandwidth smearing across the 4.9 m maximum baseline is
Δν/ν · B/λ ≈ 0.04 fringe, which is negligible.

## 6. Uncertainty

**Local covariance.** At the solution, with the analytic Jacobian *J* of
the weighted stacked residuals (*W* is inside *J*),

```
Σ_p = σ̂² (JᵀJ)⁻¹,   p = (l, m, Re C, Im C)
σ̂² = χ²/(2N − 4)          (linear loss)
σ̂² = (1.4826 MAD)²        (robust losses; insensitive to down-weighted outliers)
```

The (*l*, *m*) block, with the amplitude marginalised, is reported as
`sigma_l`, `sigma_m` and `cov_lm`. It is propagated to (Az, El) with the
analytic Jacobian

```
∂Az/∂l =  m/r²,  ∂Az/∂m = −l/r²,  ∂El/∂l = −l/(r√(1−r²)),  ∂El/∂m = −m/(r√(1−r²)),  r² = l² + m²
```

**Zenith.** When *r* < 3 × the largest (*l*, *m*) standard deviation,
azimuth is ill-defined. `ZENITH_AZ_UNCERTAIN` is set and `sigma_az_deg` is
written as NaN rather than a meaningless linearised value. Stage 4 should use
the (*l*, *m*) covariance, which stays valid there.

**Validation.** Over 60 seeded trials at per-bin SNR ≈ 7, the pulls
(error/σ) have a standard deviation of 0.8–1.0, so the uncertainties are
honest to slightly conservative. Below SNR ≈ 5 the core of the pull
distribution is still about unit width (robust std 0.98 at SNR 3.6), but
occasional fits jump to noise peaks. That is the reason for the `LOW_SNR`
threshold.

**Bootstrap.** `--bootstrap N` resamples antenna pairs with replacement
(all their samples in the bin) and re-refines from the point estimate.
Azimuth is unwrapped about the estimate, then 16/50/84 percentiles are
written (`az_boot_p16` … `el_boot_p84`). The random seed (`--seed`) is
recorded.

## 7. Signal-to-noise ratio

```
snr = |Ĉ| · √(Σw / 2) / σ_vis,      σ_vis = 1.4826 MAD of [Re r, Im r],  r = √w (V − M)
```

This is the amplitude estimate divided by its standard error. For
per-component noise σ/√w, Var(Ĉ) = 2σ²/Σw. In other words it is the
matched-filter detection significance of the point source. Because σ_vis is
taken from the residuals, it includes every unmodelled source in the sky
(confusion), not just receiver noise.

## 8. Goodness of fit

| Column | Definition |
|---|---|
| `chi2` | Σ w ‖V − M‖² |
| `reduced_chi2` | `chi2` / (2N − 4). A χ² statistic only if the weights are inverse noise variances. TART MSs have `WEIGHT = 1`, so for them it is a diagnostic in squared visibility units. |
| `rms_complex_residual` | √(Σ w ‖V − M‖² / Σ w) |
| `median_abs_residual`, `max_abs_residual` | of ‖V − M‖ |
| `sigma_vis` | robust per-component residual scale (§7) |
| `residual_tail_ratio` | √(chi2/(2N−4)) / σ_vis: about 1 for Gaussian residuals; larger means outliers or unmodelled structure |

## 9. Image-domain validation

`--image-validation` (or `--method both`/`image`) forms the real part of the
dirty image, I(l, m) = Re Σ w V e<sup>+iφ</sup> / Σ w, using the same sign
convention. The pixel size is λ/B<sub>max</sub>/6, and the image covers the
search cap. The brightest pixel inside the cap is refined by parabolic
interpolation along *l* and *m*. No CLEAN is applied.
`image_snr = peak / (1.4826 MAD)` of a hemisphere dirty image at one pixel
per resolution element. This is a dynamic-range-style measure, since its
"noise" includes all sidelobes. The disagreement is
Δθ = arccos(**s**<sub>vis</sub> · **s**<sub>img</sub>), evaluated with atan2
for precision. It is flagged `IMAGE_DISAGREEMENT` above
max(0.5 λ/B<sub>max</sub>, 3 σ<sub>position</sub>).

**Independence.** For the linear loss, minimising the profiled χ² equals
maximising the dirty-image *amplitude* ‖Σ w V e<sup>+iφ</sup>‖. The image
path differs: it uses the real part on a fixed pixel grid with interpolation,
rather than a continuous fit with a free complex amplitude and a robust
loss. It is therefore an independent implementation that catches coding and
convention errors. It is not a statistically independent measurement.

## 10. Quality flags

Several flags may apply at once, and an estimate with none is GOOD. Flags in
**bold** are blocking: they make `usable = false`.

| Flag | Condition |
|---|---|
| **`LOW_BASELINE_COUNT`** | fewer than `--min-baselines` (10) valid baselines. No fit is attempted and no direction is written. |
| `HIGH_FLAGGED_FRACTION` | more than `--max-flagged-fraction` (0.5) of the bin's cross-correlation samples were excluded |
| **`FIT_FAILED`** | the optimiser failed, or the covariance is singular or non-finite |
| **`BOUNDARY_SOLUTION`** | §4.5 |
| **`LOW_SNR`** | `snr` < `--min-snr` (5). Blocking because below SNR ≈ 5 fits can lock onto noise or sidelobe peaks (§6). |
| `HIGH_RESIDUAL` | `residual_tail_ratio` > 1.5. If the weights are calibrated (not all 1), also `reduced_chi2` > 3. |
| `IMAGE_DISAGREEMENT` | §9 |
| **`TLE_UNASSOCIATED`** | fitted direction more than `--association-radius-deg` (5°; ≤ 0 disables) from the TLE prediction: a detection of some other source. This is the post-fit association check allowed by spec §28. The fitted value is never altered. |
| `ZENITH_AZ_UNCERTAIN` | §6 |
| **`INVALID_GEOMETRY`** | solution on or below the horizon |
| `CALIBRATION_WARNING` | Stage-2 calibration issue for the station: calibration not applied, gain solutions changed during the observation, more than 25 % zero-gain antennas, or Stage-2 errors |

Excluded samples are auto-correlations, flagged samples or rows, baselines
touching an antenna flagged in the MS or a Stage-2 zero-gain antenna,
non-finite visibilities, and zero, negative or non-finite weights.

## 11. Near-field diagnostic

The plane-wave model neglects a wavefront-curvature phase of at most
π η, where η = B²/(λR). For each station it is computed with B the longest
baseline, λ at the highest channel frequency, and R the minimum
TLE-predicted range during the observation. If no TLE is available, R falls
back to the 19 100 km GNSS lower bound. The result is recorded under
`stations.<id>.near_field`.

| Station | B (m) | R (km) | η | max phase error |
|---|---|---|---|---|
| MU-UDM | 3.41 | 21 970 | 2.8 × 10⁻⁶ | 8.7 × 10⁻⁶ rad |
| ZA-Rhodes | 4.90 | 20 727 | 6.1 × 10⁻⁶ | 1.9 × 10⁻⁵ rad |

For comparison, the per-bin direction precision (σ ≈ 0.3°) corresponds to a
phase error of about 2π (B/λ) σ ≈ 0.8 rad on the longest baseline. The
curvature term is 4–5 orders of magnitude smaller, so the far-field model is
valid. The spherical-wave model was therefore **not** implemented, as the
specification requires. `far_field_valid` turns false if the error exceeds
0.01 rad.

## 12. Inputs and validation

The Stage-2 `manifest.json` is authoritative. The checks are:

- `stage == 2`;
- target `norad_id` present;
- `three_station_valid`;
- per station: `station_id` consistency, `measurement_set_status == VALID`,
  MS path, frequency and time coverage present;
- MS antenna count equal to Stage 2's;
- zero-gain antennas and calibration changes carried forward.

A dataset that is not three-station valid is refused (exit 2) unless
`--allow-invalid-stage2` is given. In that case only stations with a VALID
MS are processed, and every row and the manifest record
`stage2_three_station_valid: false`.

MS paths in the manifest are relative to the Stage-2 working directory. They
are resolved, in order:

1. relative to `--data-root`;
2. as given;
3. with the Stage-2 `run_directory` replaced by the manifest's directory;
4. against each ancestor of the manifest directory.

`--station` accepts a `station_id` (`mauritius`) or a `tart_id` (`mu-udm`).

## 13. Time binning

Bins of `--time-bin` seconds (default 1) start at each station's first
timestamp. Empty bins are omitted, so no timestamp is ever fabricated. The
output time is the weight-averaged time of the samples actually used. TART
cadence is about 1.025 s, so 1 s bins usually hold one dump.

## 14. Outputs

```text
results/stage3/20260625/NORAD_40534/
├── stage3_manifest.json
├── directions.csv        # quality_flags ';'-separated
├── directions.parquet    # quality_flags as a list (needs pyarrow)
├── diagnostics/<station>_{track,sky,vis_image_offset,fit_quality}.png
└── logs/stage3.log
```

Table columns (schema version 1, one row per station and non-empty bin):

- **Identity:** `time_utc`, `time_mjd` (UTC MJD, days), `station_id`,
  `tart_id`, `norad_id`, `primary_method`, `n_times`, `n_samples`.
- **Direction:** `l`, `m`, `sigma_l`, `sigma_m`, `cov_lm`, `az_obs_deg`,
  `el_obs_deg`, `sigma_az_deg`, `sigma_el_deg`, `cov_az_el`.
- **Amplitude:** `complex_amp_real`, `complex_amp_imag`,
  `apparent_amplitude`.
- **Fit quality:** `chi2`, `reduced_chi2`, `rms_complex_residual`,
  `median_abs_residual`, `max_abs_residual`, `sigma_vis`,
  `residual_tail_ratio`, `snr`.
- **Data selection:** `n_baselines_total`, `n_baselines_used`,
  `flagged_fraction`.
- **Image cross-check:** `az_image_deg`, `el_image_deg`, `image_peak_value`,
  `image_snr`, `vis_image_separation_deg`.
- **Bootstrap:** `az_boot_p16|p50|p84`, `el_boot_p16|p50|p84`.
- **Status:** `fit_success`, `optimizer_status`, `usable`, `quality_flags`,
  `stage2_three_station_valid`.
- **Reference predictions** (not measurements): `az_tle_deg`, `el_tle_deg`,
  `range_tle_km`, `tle_separation_deg`.

`stage3_manifest.json` contains:

- `schema_version`, `stage: 3`;
- the input manifest path and SHA-256;
- `norad_id` and the satellite and TLE;
- per-station summaries: calibration, UVW check, near-field result, counts,
  fit-success fraction, median SNR, median visibility–image offset, flag
  counts;
- every estimator option, including robust loss, search radius, grid,
  uncertainty and SNR definitions, bootstrap count and seed;
- time binning and quality thresholds;
- software versions (tart_gnss, Python, NumPy, SciPy, pandas, pyarrow,
  matplotlib, python-casacore) with the git commit and dirty flag;
- the command line, a UTC timestamp, and the output paths.

The figures are generated from the saved table, not from in-memory objects.

**Exit status:** `0` means at least one usable estimate; `3` means it
completed with none; `2` means invalid input or configuration (manifest,
station, missing MS, convention mismatch); `1` means an unexpected failure.

## 15. Commands

```bash
pip install -e ".[stage3]"

# all stations, visibility estimator + image cross-check
tart-localise --manifest data/observations/20260625/NORAD_40534/manifest.json \
    --all-stations --method visibility --image-validation --time-bin 1 \
    --output-dir results/stage3

# one station by TART id, 10 s bins, soft-L1 loss, 200 bootstrap resamples
tart-localise --manifest data/observations/20260625/NORAD_40534/manifest.json \
    --station mu-udm --time-bin 10 --robust-loss soft_l1 --bootstrap 200 \
    --output-dir results/stage3

# manifest written in another working directory
tart-localise --manifest /path/to/manifest.json --all-stations --data-root /path/to/stage2/workdir

# synthetic end-to-end demonstration (writes docs/examples/stage3/)
python scripts/stage3_synthetic_demo.py
```

## 16. Synthetic end-to-end demonstration

[`docs/examples/stage3/`](examples/stage3) holds the demonstration
([`scripts/stage3_synthetic_demo.py`](../scripts/stage3_synthetic_demo.py)).
It uses three synthetic 24-antenna arrays at the real station coordinates,
observing for 20 s on 2026-06-25 from 09:14 UTC. The source was injected
3.0–3.7° away from the real TLE track of GPS PRN 26 (more than 20σ), with
per-component noise 0.5, |C| = 0.36 and 5 % random flags. Values below are for the middle
bin; the full table is in [`demo_report.txt`](examples/stage3/demo_report.txt).

| Station | Injected Az/El | Recovered Az/El | Error | σ<sub>Az</sub>/σ<sub>El</sub> | Image Az/El | Vis–img | TLE sep. |
|---|---|---|---|---|---|---|---|
| Mauritius | 217.563 / 45.014 | 217.689 / 45.189 | 0.196° | 0.144 / 0.129° | 217.438 / 45.462 | 0.33° | 3.39° |
| Madagascar | 199.878 / 50.191 | 199.710 / 50.180 | 0.108° | 0.159 / 0.123° | 199.599 / 50.339 | 0.17° | 3.70° |
| South Africa | 181.449 / 68.184 | 181.299 / 68.216 | 0.064° | 0.307 / 0.189° | 180.804 / 67.883 | 0.38° | 3.21° |

Over all 60 bins, the median error is 0.11–0.17° (maximum 0.38°) and the
median σ<sub>El</sub> is 0.13–0.17°. The recovered direction stays
3.0–3.7° from the TLE, so it follows the injected source, not the search
centre. Every bin is GOOD.

## 17. Real-data run: 2026-06-25, GPS PRN 26 (NORAD 40534)

Input: the valid Stage-2 dataset (three 24-antenna MSs, 645 s real overlap).
The Measurement Sets are not in git. The run was
`tart-localise --all-stations --image-validation --time-bin 10 --bootstrap 50`
in the Stage-2 working directory, with exit status 0. TLE: epoch 103 days
after the observation (reference only).

| Station | Bins | Usable | Usable interval (UTC) | Median TLE sep. | Median ΔAz / ΔEl vs TLE | σ<sub>Az</sub>/σ<sub>El</sub> (Jacobian) | Bootstrap half-width Az/El | SNR | Vis–img |
|---|---|---|---|---|---|---|---|---|---|
| ZA-Rhodes (15/24 antennas) | 68 | **68** | 09:13:13–09:24:21 | 0.96° | +0.10° / −0.92° | 0.17° / 0.10° | 0.97° / 0.68° | 12.2 | 0.31° |
| MU-UDM | 68 | **37** | 09:18:36–09:24:34 | 1.42° | −1.96° / +0.12° | 0.17° / 0.15° | 0.55° / 0.48° | 9.5 | 0.28° |
| MG-Tana | 68 | **2** | 09:23:54–09:24:43 | 2.43° | −3.59° / +0.06° | 0.19° / 0.15° | 0.83° / 0.45° | 9.6 | 0.26° |

What the flags show:

- **ZA-Rhodes** tracks PRN 26 for the whole interval, 1° from the stale TLE
  track. Every bin carries `HIGH_FLAGGED_FRACTION` (62 % of samples lost to the
  9 zero-gain antennas) and `CALIBRATION_WARNING`; 8 carry `HIGH_RESIDUAL`.
- **MU-UDM**: in the first half (09:13–09:18) a brighter source about 10–13°
  from the TLE direction pulls the fit to the cap edge (31 bins
  `BOUNDARY_SOLUTION` and `TLE_UNASSOCIATED`). From 09:18:36 the fit locks
  onto the PRN 26 track.
- **MG-Tana**: the fit is dominated by non-moving sources: about
  Az 193.6°, El 47.6° at the start and end of the run, and about Az 197°,
  El 37–38° (at the cap edge) in between. They stay put for 11 minutes
  while PRN 26's predicted track moves 7° in azimuth, so they are **not**
  PRN 26. It is presumably a near-geostationary or terrestrial
  emitter. 60 bins are flagged `TLE_UNASSOCIATED` and 55
  `BOUNDARY_SOLUTION`. PRN 26 itself is only seen at the very end.
- **1 s bins** (the default): per-bin SNR is 3.0–3.9 (median), so 1979 of
  1980 bins are `LOW_SNR`. Use `--time-bin 10` for TART. The satellite moves
  about 0.3° in 10 s, well within a 3° beam, and the SNR rises by about √10.

**Uncertainty on real data.** The Jacobian σ was checked against the
data themselves. 10 s bins do not overlap, so their scatter is an
independent estimate of the random error:

| Station (usable bins) | Bin-to-bin scatter (successive differences / √2), El / Az·cos El | Residual about a cubic track, El / Az·cos El | Lag-1 autocorrelation of El residuals | Jacobian σ, El / Az·cos El | Bootstrap half-width, El / Az·cos El |
|---|---|---|---|---|---|
| ZA-Rhodes (68) | 0.038° / 0.022° | 0.177° / 0.063° | 0.96 | 0.104° / 0.062° | 0.68° / 0.37° |
| MU-UDM (37) | 0.048° / 0.041° | 0.071° / 0.066° | 0.61 | 0.153° / 0.118° | 0.48° / 0.39° |

What this shows:

- **The per-bin Jacobian σ is honest to conservative** for the random
  (bin-to-bin) error.
- **The bootstrap overstates it by about 10×.** Resampling the baselines of
  a sparse 15–24-element array changes its sidelobe pattern, so on a crowded
  sky the bootstrap measures sensitivity to confusion, not measurement noise.
- **The dominant real-data error is systematic.** It is slowly varying and
  correlated over minutes: there are 0.1–0.2° wanders about a smooth track
  (high autocorrelation), and the 1–2° offsets from the TLE track are
  consistent with the 103-day TLE propagation, residual calibration errors
  and sidelobe confusion. Stage 3 cannot separate these.

Stage 4 should therefore not treat successive bins as independent, and
should model a per-station systematic term rather than rely on either σ
alone.

The near-field metric is η = 2.8 × 10⁻⁶ (MU), 2.5 × 10⁻⁶ (MG) and
6.1 × 10⁻⁶ (ZA), so the plane-wave model holds (§11).

## 18. Limitations

- **Single point source.** The model has one source. TART sees 10–15 GNSS
  satellites at once, and their sidelobes bias the fit and inflate σ_vis.
  The search cap does the source association. A multi-source or
  subtract-and-refit model is future work.
- **Sidelobe capture.** TART's sparse 24-element array has high sidelobes.
  A brighter source outside the cap can pull the solution to the cap edge
  (flagged). A source beyond the guard band can go undetected, leaving a
  sidelobe solution inside.
- **Calibration** is taken as applied by `tart2ms` (Stage 2: *applied,
  unverified*). Gain-solution changes and zero-gain antennas are flagged, not
  corrected. Residual antenna phase errors bias directions at a level not
  captured by σ.
- **Weights** in TART MSs are all 1, so `reduced_chi2` is not a calibrated χ²
  and σ is estimated from the residuals.
- **Uncertainties** are local (Gauss–Newton) per-bin random errors. They are
  honest on synthetic data at SNR ≳ 5, and match the bin-to-bin scatter of
  real data. They do **not** include the slowly varying systematics
  (calibration, confusion) that dominate real TART errors and are correlated
  across bins (§17). On a crowded sky the baseline bootstrap measures
  confusion sensitivity and overstates the random error.
- **Image estimator** is not statistically independent of the visibility
  estimator (§9).
- **TLE reference.** The TLE in the 2026-06-25 example has an epoch 103 days
  after the observation. Its track is a rough reference, which matters only
  for the search centre and the plots.
- **Time bins** are per station. Nothing is interpolated or aligned across
  stations (Stage 4).
