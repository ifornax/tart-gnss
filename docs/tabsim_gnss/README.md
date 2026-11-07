# GNSS tab-sim validation: reproducible Ghana experiment

[Technical LaTeX report](../report/tabsim_gnss_validation.tex) ·
[Orbit-blind fitting script](../../scripts/tabsim_gnss_fit.py) ·
[MS and orbit diagnostics](../../scripts/tabsim_gnss_diagnostics.py)

This is a research validation of synthetic visibilities, not an accuracy claim
for real TART data. Existing Stage 1–7 processing is unchanged.

## Environment

Run in the local tart-gnss checkout, with the tabascal environment:

~~~bash
conda activate tabascal
export REPO="$HOME/GitHub/TART/sat_position"
export TART_CONFIG="$REPO/configs/tart_tabsim"
export TART_SIM="$TART_CONFIG/simulated"
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

## Run simulation (when not already generated)

The three site-specific YAMLs live outside the repository at TART_CONFIG.
They specify NORAD 39166, GNSS L1 1575.42 MHz, 24 antennas,
120 integrations of 1.024 s, 33 samples per integration, and a
1 MHz channel. Ensure output paths do not collide with existing
runs if overwrite=false.

~~~bash
sim-vis -c "$TART_CONFIG/tab_sim_ghana.yaml"
sim-vis -c "$TART_CONFIG/tab_sim_mg-tana.yaml"
sim-vis -c "$TART_CONFIG/tab_sim_za-rhodes.yaml"
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
claiming identical orbital truth. Never commit large MS/Zarr outputs,
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
