# GNSS tab-sim validation: reproducible Ghana experiment

[Technical LaTeX report](../report/tabsim_gnss_validation.tex) ·
[Orbit-blind fitting script](../../scripts/tabsim_gnss_fit.py) ·
[MS and orbit diagnostics](../../scripts/tabsim_gnss_diagnostics.py) ·
[Monte Carlo uncertainty](../../scripts/tabsim_ghana_uncertainty.py) ·
[Frame/time audit](../../scripts/tabsim_ghana_frame_time_audit.py) ·
[UTC/UT1 check](../../scripts/tabsim_ghana_time_scale_check.py) ·
[Truth re-scoring](../../scripts/tabsim_ghana_rescore_truth.py)

This is a research validation of synthetic visibilities, not an accuracy claim
for real TART data. Existing Stage 1–7 processing is unchanged.

## Environment

Run from the root of a tart-gnss checkout, with the tabascal environment.
The simulation inputs are committed in
[`configs/tart_tabsim/`](../../configs/tart_tabsim) (see
[Simulation configuration files](#simulation-configuration-files)). The
simulated Measurement Sets and Zarr stores are large and stay local. The
YAML files write them to `outputs/tabsim/<station>/` (git-ignored). The
original runs reported here were written under the author's `sat_position`
working copy; if you already have runs elsewhere, point `TART_SIM` there:

~~~bash
conda activate tabascal
export REPO="$(git rev-parse --show-toplevel)"                    # this repository (run inside the checkout)
export TART_CONFIG="$REPO/configs/tart_tabsim"                    # committed inputs
export TART_SIM="$REPO/outputs/tabsim"                            # local outputs (override if your runs live elsewhere)
cd "$REPO"
python -m pip show tabsim
python -m pip check
python -c "import astropy, skyfield, scipy, casacore.tables, numpy; print('OK')"
~~~

The successful simulator was tabsim 0.0.13, from
https://github.com/chrisfinlay/tab-sim. The original local
configuration generator is
$HOME/GitHub/SBI_AI/tart_tabsim_three_station/make_configs.py,
which derives station coordinates from actual MS antenna tables.
Do not substitute invented station layouts or SEFD values.

## Simulation configuration files

Copied from the local simulation working tree
(`sat_position/configs/tart_tabsim/`, 2026-11-17). The top-level
`tab_sim_*.yaml` inputs differ from the original runs only in
`itrf_path` and `output.path`, which are now relative to the repository
root so that they work on any machine. The exact as-run record is
`simulated/<station>/<run>/input_data/sim_config.yaml`:

| File | Content |
|---|---|
| `tab_sim_{ghana,mg-tana,za-rhodes}.yaml` | tab-sim input for each station: site latitude/longitude/elevation, ITRF antenna file, NORAD 39166 phase centre (RA/Dec), start 2026-10-29T01:10:00 UTC, 120 × 1.024 s, 32 sub-integrations, 1575.42 MHz, 1 MHz channel, SEFD 10 000 Jy, unit gains, seeds |
| `{ghana,mg-tana,za-rhodes}.itrf.txt` | 24 antenna positions per station (ITRF/ECEF, metres), derived from the real MS antenna tables |
| `simulated/<station>/<run>/input_data/sim_config.yaml` | the full configuration tab-sim resolved for that run, with every default filled in |
| `simulated/<station>/<run>/input_data/norad_ids.yaml`, `*.itrf.txt` | the NORAD list and antenna file as used by that run |
| `simulated/<station>/<run>/input_data/used_orbits.json` | the exact OMM record tab-sim propagated (NAVSTAR 68, NORAD 39166, CelesTrak epoch 2026-10-30T07:09:07.795296 UTC). The scripts read it via `--orbits` to score fits. The three runs carry identical orbital elements and differ only in the `FETCHED_AT` stamp (absent for Ghana) |
| `simulated/<station>/<run>/log_sim_*.txt` | tab-sim run log: OMM record used (SatChecker/CelesTrak epoch 2026-10-30T07:09:07.795 UTC, age 1.25 d), source altitude, write times |

The SEFD, `dish_d`, `chan_width` and `power_scale` values are exploratory
inputs, not measured TART properties. The run directories hold only these
text files and `used_orbits.json` in git; the `.ms`, `.zarr`, plots and
spectral-model files are git-ignored.

**Paths.** tab-sim uses `itrf_path` and `output.path` exactly as written:
it does not expand `~` or environment variables, and it resolves relative
paths against the current directory. Run `sim-vis` from the repository
root. Only the as-run records (`sim_config.yaml`, `log_sim_*.txt`) keep
the original machine's absolute paths, as provenance; nothing reads them.

## Run simulation (when not already generated)

The three site-specific YAMLs are in `configs/tart_tabsim/`. Run them
from the repository root; each writes to `outputs/tabsim/<station>/`. They specify NORAD 39166, GNSS L1 1575.42 MHz, 24 antennas,
120 integrations of 1.024 s, 33 samples per integration, and a
1 MHz channel. Ensure output paths do not collide with existing
runs if overwrite=false.

~~~bash
cd "$REPO"
sim-vis -c configs/tart_tabsim/tab_sim_ghana.yaml
sim-vis -c configs/tart_tabsim/tab_sim_mg-tana.yaml
sim-vis -c configs/tart_tabsim/tab_sim_za-rhodes.yaml
~~~

Confirm the output directory has 1SAT, not 0SAT. The failed old
simulation contained noise only. The successful simulation used
the saved OMM epoch 2026-10-30T07:09:07.795296 UTC and a
one-satellite spectral emission model. All three observed 1SAT
files have nonzero RFI_MODEL_DATA.

Select the correct Ghana MS explicitly if there are multiple runs:

~~~bash
export GHANA_MS="$(find "$TART_SIM/ghana" -type d -name '*_1SAT_*.ms' -print | sort | tail -1)"
test -d "$GHANA_MS" || { echo "Ghana MS missing"; exit 1; }
export GHANA_RUN="$(dirname "$GHANA_MS")"
export GHANA_ORBITS="$GHANA_RUN/input_data/used_orbits.json"
echo "$GHANA_MS"
~~~

The same orbit record is committed, so `GHANA_ORBITS` can also point at
`$TART_CONFIG/simulated/ghana/*/input_data/used_orbits.json`.

## Inspect MS and validate forward-model phases

~~~bash
python scripts/tabsim_gnss_diagnostics.py inspect --ms "$GHANA_MS"
python scripts/tabsim_gnss_diagnostics.py phase \
    --ms "$GHANA_MS" --orbits "$GHANA_ORBITS" --samples 33
~~~

The second command deliberately uses **orbital truth** for forward-model
validation, not fitting. Results observed on Ghana:

| Phase model | Coherence | Circular RMS |
|---|---:|---:|
| Finite range plus antenna-w | ~1.000000 | 0.0479 degrees |
| Finite range without w | 0.037762 | 103.4545 degrees |
| Conjugated corrected model | 0.826855 | 35.0277 degrees |

## Fit direction without orbit truth

~~~bash
mkdir -p results/tabsim_validation
python scripts/tabsim_gnss_fit.py "$GHANA_MS" \
    --snapshot 0 --range-km 21500 --samples 33 \
    --search-deg 15 --seed 123 --column DATA \
    --output results/tabsim_validation/ghana_direction_fit.json
~~~

The orbit-blind fit uses MS FIELD pointing to set a bounded local
search, and eliminates a single complex nuisance amplitude analytically.
It does not read the saved OMM. In the observed idealised simulation:

| Parameter | Ghana result |
|---|---:|
| Azimuth | 100.948434157 deg |
| Elevation | 51.728550893 deg |
| Normalised squared residual | 5.750019e-12 |
| Phase coherence | 0.999999999997 |
| Local optimiser | success |
| Global optimiser | **failure status** |

An independent OMM prediction at 2026-10-29T01:10:00 UTC:
Az = 100.950025 deg, El = 51.728246 deg and range = 21498.801 km.
The corresponding internal angle agreement is roughly 3.7 arcsec,
not a measured astronomical localisation accuracy.

## Nominal-range sensitivity

~~~bash
for range in 18000 20000 21500 23000 26000; do
  python scripts/tabsim_gnss_fit.py "$GHANA_MS" \
    --range-km "$range" --snapshot 0 --samples 33 \
    --search-deg 15 --seed 123 \
    --output "results/tabsim_validation/ghana_fit_"$range"km.json"
done
~~~

The five observed fits show at most ~0.00032 arcsec variation
in fitted Az and ~0.00009 arcsec in El; their objective values
are nearly identical. Single-station range is weakly identified
and the differences may be numerically insignificant.

## Compile LaTeX

~~~bash
cd "$REPO/docs/report"
pdflatex -interaction=nonstopmode -halt-on-error tabsim_gnss_validation.tex
pdflatex -interaction=nonstopmode -halt-on-error tabsim_gnss_validation.tex
~~~

Requires standard LaTeX packages: geometry, amsmath, amssymb, booktabs,
graphicx, hyperref, listings, xcolor, longtable.

## Provenance and open work

Orbit records are saved per run as input_data/used_orbits.json;
archive records from different runs must be compared in full before
claiming identical orbital truth. Only the YAML/text inputs and run
records in `configs/tart_tabsim/` (including `used_orbits.json`) are committed. Never commit large MS/Zarr outputs,
credentials, or machine-specific files. We have not independently
measured the instrument SEFD or calibrated angular uncertainty.

**Next validation:** multiple seeds and global minima; realistic
noise/power and gain errors; empirical covariance and 68% interval
coverage; direction fitting at Antananarivo and Rhodes; then compare
three-site Bayesian triangulation with existing Stage 4.

**Decision log (2026-11-07):** adopt tab-sim as the controlled forward
simulator; use phase-centre-corrected spherical-wave delays; keep
two-stage direction → 3D inference as a working hypothesis pending
uncertainty calibration. This report does not modify operational stages.


## Ghana directional-uncertainty Monte Carlo (next experiment)

The new [uncertainty script](../../scripts/tabsim_ghana_uncertainty.py) generates
independent synthetic noise perturbations of the existing satellite
`RFI_MODEL_DATA`, fits each perturbed snapshot with the same orbit-blind
finite-range + \(w\)-corrected forward model, and uses saved-orbit truth
**only after the fitting loop** to assess accuracy.

The per-visibility SNR is defined as

\[
{\rm SNR}_{\rm vis} =
\frac{\sqrt{\langle |S_b|^2\rangle_b}}
{\sqrt{\langle |N_b|^2\rangle_b}}
=\frac{\sqrt{\langle |S_b|^2\rangle_b}}{\sqrt{2}\,\sigma},
\]

where independent noise is generated as
\(N_b=\sigma(Z_{1b}+iZ_{2b})\), with independent standard normal
\(Z_{1b},Z_{2b}\). The value of \(\sigma\) is estimated from the
full MS `NOISE_DATA` (all unflagged samples), while the template
\(S_b\) comes from the selected `RFI_MODEL_DATA` snapshot.
This is not the coherent array detection SNR.

### Initial smoke test

Run from repository root after defining `GHANA_MS` and
`GHANA_ORBITS` above:

~~~bash
python scripts/tabsim_ghana_uncertainty.py \
  --ms "$GHANA_MS" \
  --orbits "$GHANA_ORBITS" \
  --snrs 3,10,30 \
  --reps 12 \
  --seed 20261107 \
  --samples 33 \
  --range-km 21500 \
  --search-deg 3 \
  --coarse-step-deg 1 \
  --output-dir results/tabsim_validation/ghana_uncertainty_smoke
~~~

The small trial count is a programmatic smoke test **only**, not a
stable covariance estimate. The fitter uses a bounded pointing-informed
grid search plus local Nelder-Mead; monitor `local_success` and
`on_edge` flags before trusting results.

### Full experiment

~~~bash
python scripts/tabsim_ghana_uncertainty.py \
  --ms "$GHANA_MS" \
  --orbits "$GHANA_ORBITS" \
  --snrs 3,5,10,20,50,100 \
  --reps 200 \
  --seed 20261107 \
  --samples 33 \
  --range-km 21500 \
  --search-deg 3 \
  --coarse-step-deg 1 \
  --output-dir results/tabsim_validation/ghana_uncertainty
~~~

Outputs:

- `ghana_localisation_error_vs_snr.png`: median, 68th and 95th
  percentile sky-plane angular error as a function of input SNR.
- `ghana_direction_trials.csv`: individual fits and ENU tangent-plane
  errors, retaining optimisation success/boundary flags.
- `ghana_direction_snr_summary.csv`: per-level RMSE, quantiles,
  bias and independent held-out 68% ellipse coverage.
- `ghana_direction_covariance.json`: covariance matrices in
  **East/North arcsec squared**, separated by injected SNR.
- `run_metadata.json`: all experiment parameters and RNG seed.

The East/North errors are
\(e_E=3600\cos(El_\star)\operatorname{wrap}(Az-Az_\star)\) and
\(e_N=3600(El-El_\star)\), in arcseconds. Each SNR group is split:
first 50% to estimate mean bias and covariance, remaining 50%
to check joint-Gaussian 68% ellipse coverage
(\(D^2\le\chi^2_{2,0.68}\approx2.279\)).

**Triangulation caution:** \(\Sigma_{EN}\) here is empirical
estimator scatter **conditional on one synthetic source direction
and one injected noise model**. It is not automatically the
posterior covariance of an individual direction estimate.
Check calibration and bias, then convert to
\(\Sigma_{Az,El} = D^{-1}\Sigma_{EN}(D^{-1})^\mathsf{T}\),
\(D=\operatorname{diag}(3600\cos El_\star,3600)\)
for degree-squared covariance if Stage 4 accepts degrees.
Avoid treating held-out coverage estimated using an uncorrected
systematic bias as the same quantity as bias-corrected coverage.

The full experiment (200 trials per SNR, 1,200 fits) has since been
run on the local Ghana MS. Its results, re-scored against both time
references, are in
[the re-scoring section](#rescore-the-existing-monte-carlo-runs-after-the-utcut1-diagnosis).


## Diagnosing the persistent 3.713-arcsec noise-free offset

The earlier Ghana audit found a deterministic noise-free fit residual of
(E, N) = (-3.548, +1.096) arcsec. It remained stable when changing
the fixed nominal range 18,000--26,000 km and grid spacing 1.0--0.5 deg.
The unsuccessful Nelder-Mead status flags at finite SNR did not
correspond to obviously worse fits, but termination details still need
to be recorded.

To distinguish a coordinate-frame mismatch from the moving-source
approximation, run the new frame and time audit:

~~~bash
python scripts/tabsim_ghana_frame_time_audit.py \
  --ms "$GHANA_MS" \
  --orbits "$GHANA_ORBITS" \
  --samples 33 --range-km 21500 \
  --search-deg 3 --coarse-step-deg 0.5 \
  --output results/tabsim_validation/ghana_audit/frame_time_bias.json
~~~

The report compares (a) topocentric Skyfield Az/El,
(b) orbital Az/El computed with exactly the estimator's GCRS
antenna-centroid ENU basis, (c) time-mean orbital direction over
the integration, and (d) the noise-free fitted direction.
It also evaluates matched near-field visibilities under a moving-orbit
model and stationary-direction models.

Interpretation guide: if Skyfield and GCRS/ENU truth differ by
about 3.7 arcsec, a reference-frame or observer-location mismatch
is implicated. If the frame values agree but the mean orbital
direction shifts toward the fit, motion/integration conventions
are implicated. If neither matches, compare time-dependent
station antenna-w coordinates, beam weighting, and orbital/GCRS
propagation conventions before applying a bias correction.

The script is a diagnostic: it does not modify or recalibrate
the Monte Carlo covariance results. Running it, together with the
controlled UTC/UT1 test below
(`scripts/tabsim_ghana_time_scale_check.py`), traced the offset to the
simulator's time-scale convention, not to a frame mismatch:

~~~bash
python scripts/tabsim_ghana_time_scale_check.py \
  --ms "$GHANA_MS" --orbits "$GHANA_ORBITS" --samples 33 \
  --output results/tabsim_validation/ghana_audit/time_scale_check.json
~~~


## Rescore the existing Monte Carlo runs after the UTC/UT1 diagnosis

A controlled Skyfield phase test using identical antenna geometry showed:

- physical UTC: profile residual \(J=3.82488077454\times10^{-7}\);
- simulator-compatible UT1: \(J=1.20662429051\times10^{-12}\);
- simulated minus physical UTC reference: East \(-3.54816\) arcsec,
  North \(+1.09650\) arcsec.

This exactly explains the noise-free \(\approx 3.713\)-arcsec
reference discrepancy. The fit itself recovers the direction
encoded by the simulated phases. For this run, `tabsim` propagates
numeric observation Julian dates through `ts.ut1_jd`, although the
observation calendar date is recorded as UTC. Document both references;
do not silently redefine UTC as UT1.

Rescore the original 1,200 trial records **without rerunning fits**:

~~~bash
python scripts/tabsim_ghana_rescore_truth.py \
  --trials results/tabsim_validation/ghana_uncertainty/ghana_direction_trials.csv \
  --time-check results/tabsim_validation/ghana_audit/time_scale_check.json \
  --output-dir results/tabsim_validation/ghana_truth_corrected
~~~

If `time_scale_check.json` is absent, explicit reference values from
the verified diagnostic can be supplied via `--sim-az`, `--sim-el`,
`--utc-az`, and `--utc-el`. This writes parallel
`simulator_ut1` and `physical_utc` subdirectories containing
trial CSV, per-SNR summaries, covariance JSON and error plots.

Observed rescoring of the supplied 200-trial-per-SNR CSV:

| SNR | UTC median (arcsec) | Simulator median (arcsec) | Simulator raw 68% held-out coverage |
|---:|---:|---:|---:|
| 3 | 110.72 | 111.45 | 55% |
| 5 | 63.62 | 63.97 | 71% |
| 10 | 30.54 | 30.43 | 74% |
| 20 | 16.05 | 16.48 | 54% |
| 50 | 7.05 | 6.55 | 71% |
| 100 | 4.90 | 3.35 | 67% |

At SNR 100, training-set East/North bias reduces from
approximately (-3.57,+0.88) to (-0.03,-0.22) arcsec. Corrected
covariances are numerically almost unchanged because they are centered
on the training errors; the angular coordinate scale changes slightly
with the chosen reference elevation.

**Upstream simulator fix is not yet applied.** In the installed
`tabsim` code, investigate the observation times in
`tabsim/dask/observation.py` and transforms in
`tabsim/dask/coordinates.py`, together with
`tabsim/tle.py:get_satellite_positions`. A candidate resolution is
to convert UTC timestamp values to *actual UT1 JD* using an authoritative
time-scale conversion before passing them to the Skyfield `ts.ut1_jd`
functions, or to use Skyfield UTC datetimes consistently. This must
apply identically to satellite and antenna propagation, not only to
the truth-evaluation script. Add a regression test comparing orbit
Az/El and baseline phases at a known UTC time across the simulator and
an independent implementation. Do not patch a simulator installation
or regenerate production MS files before those tests pass.
