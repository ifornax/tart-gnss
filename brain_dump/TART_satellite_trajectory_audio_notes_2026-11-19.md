# TART Satellite Trajectory — Audio Transcription and Project Notes

**Source:** `Sbi_tart.m4a` (approximately 2 min 30 s)  
**Context:** SBI workshop discussion, November 2026  
**Status:** Working project note; not an implementation decision or validation result.

## 1. Edited transcription

> In this project, what we are trying to do is to understand, or to characterise, uncertainty in a satellite's position.
>
> What we have done so far is, using TABASCAL simulation, we simulated the location of one satellite for a TART station.
>
> The next step that we did was to write down the visibility as measured by antennas p and q, and carry out the forward model to simulate—to calculate—the visibility measured for each baseline.
>
> The plan is to find out the uncertainty from the simulated TABASCAL visibility and the forward model.
>
> One caveat here is that, in reality, we use satellite TLEs to calibrate the antennas. We may have multiple satellites in the field of view that are used to calibrate the antennas.
>
> But what we need to do with real data is to find one satellite that is common between two stations, and then exclude it from the calibration. We use all the other satellites in the field of view to calibrate a station—let's say Mauritius—excluding the common satellite between Mauritius and, for example, Ghana.
>
> Afterwards, what we will try to do is characterise the position of that common satellite.

**Editorial note:** Background music, incidental speech, pauses and repeated words were omitted. The recording refers to 'satellite TLEs' in the calibration discussion; this should be checked against the actual calibration implementation rather than assumed as an established algorithm.

## 2. Project interpretation and implementation notes

### Scientific goal

Estimate the position or trajectory of a satellite observed by geographically separated TART stations, and quantify positional uncertainty from interferometric visibility measurements.

### Work described as completed in the recording

1. Generate simulated visibilities for a single satellite at a TART station using TABASCAL.
2. Formulate a baseline-dependent visibility forward model for antenna pairs `(p, q)`.

**Validation not yet established by this recording:** Numerical agreement between TABASCAL and the independent forward model, uncertainty calibration, and recovery of the satellite's true position.

### Proposed real-data calibration strategy

1. Identify a satellite simultaneously visible to two TART stations (Mauritius and Ghana are the example in the recording).
2. Exclude that target satellite from the calibration model to avoid target-information leakage.
3. Calibrate each station using the remaining suitable satellites in its field of view.
4. Apply the resulting calibration to target-satellite visibility data.
5. Infer the target satellite's position and its uncertainty, using an appropriate likelihood, prior and forward model.

### Technical cautions and unresolved decisions

- Clarify exactly how satellite TLEs enter calibration, and which station-dependent instrumental parameters are solved.
- A common satellite in two stations' fields of view does not automatically imply identical visibility time sampling or frequency coverage; establish overlap explicitly.
- Confirm whether the initial inference is independent per station or jointly uses multiple stations, and whether the target is direction, range, full 3D position, or trajectory parameters.
- Define the noise/covariance model and uncertainty method. Comparing two deterministic visibility predictions alone does not yield a calibrated positional uncertainty.
- Ensure the target is genuinely excluded from calibration, including any sky or satellite model used in gain solving.
- Verify TLE propagation accuracy and epoch validity when using it as a prior or simulation truth.

## 3. Suggested project tracking entry

**Issue:** Target-excluded calibration and satellite-position uncertainty estimation

**Acceptance criteria (proposed, not yet approved):**

- Demonstrate forward-model agreement with TABASCAL simulated visibilities at fixed satellite truth and documented tolerances.
- Show calibration using non-target satellites while holding out the common target satellite.
- Define a likelihood and inference procedure that recovers satellite parameters and credible intervals on simulated data.
- Check bias, interval coverage, and sensitivity to calibration errors and TLE prior width.
- Document data provenance, station configuration, assumptions, results, and open questions.

## 4. Source and provenance

Transcribed from the audio provided in this conversation. This document is a draft for inclusion in the TART satellite trajectory project's documentation; no remote repository has been changed.
