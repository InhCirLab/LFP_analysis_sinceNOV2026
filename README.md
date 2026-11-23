# LFP_analysis_sinceNOV2026

MATLAB pipeline for hippocampal LFP analysis (theta/gamma during locomotion, 1–120 Hz power spectra per layer, and sharp-wave ripple detection/characterization) from silicon-probe recordings.

Maintained by the **L234 lab, CNC**.

> **Important:** this repo does **not** contain the preprocessing step. First run the preprocessing from [HippoCookBook](https://github.com/valegarman/HippoCookBook) (Manu's code). The scripts here start from its output.

---

## Table of contents

1. [What this repo does](#what-this-repo-does)
2. [Repository structure](#repository-structure)
3. [Requirements](#requirements)
4. [Expected input data](#expected-input-data)
5. [How to run the main pipeline](#how-to-run-the-main-pipeline)
6. [Outputs](#outputs)
7. [Key parameters](#key-parameters)
8. [Optional analyses](#optional-analyses)
9. [Known issues / TODO](#known-issues--todo)
10. [Contact](#contact)

---

## What this repo does

`LFPpipeline.m` is the main script. It is organised in parts that must be run in this order:

| Step | Part | What it does |
|------|------|--------------|
| 1 | **P1 – Pass 1** | Theta detection + locomotion (run) epochs. Only generates `thetaRun.newintervals`, which the CNN needs to exclude locomotion. No Excel output. |
| 2 | **P3 – Ripples** | CNN-based ripple detection and classification per group, ripple measures, per-subject ID, hierarchical robust outlier removal, histograms. Saves the **new** `<rec>.ripples.events.mat`. |
| 3 | **P1R – Refresh** | Re-runs theta detection, now reading the new ripples, so `intervals_nonTheta_nonRipples` is up to date. Rewrites `.thetaEpochs6to12.states.mat`, saves Excel files and figures. |
| 4 | **P2 – Power spectrum** | 1–120 Hz power spectra per hippocampal layer (`ps1to120_Paria`), for delta/theta/beta/low-, medium- and high-gamma, in several conditions. |

**Why this order?** `intervals_nonTheta_nonRipples` is computed from whatever `<rec>.ripples.events.mat` exists at that moment. Running ripples (P3) *before* the refresh (P1R) and the power spectrum (P2) guarantees that non-run and non-theta/non-ripple epochs use the new ripples, not the old ones.

The script runs **serially** (no parallel pool) to make debugging easier. See the header of `LFPpipeline.m` for how to switch back to `parfor`.

---

## Repository structure

```
LFPpipeline.m                      Main pipeline (P1 → P3 → P1R → P2)
README.md

mandatory_functions/               Functions required by LFPpipeline.m
  runThetaStage.m                  Runs the theta/run stage over all groups and sessions
  processRecording.m               Processes ONE recording (theta/run/non-run)
  detectThetaEpochs.m              Theta epoch detection
  getWavelet.m                     Wavelet transform helper
  saveStates.m                     Saves the .states.mat files
  getAnimalID.m                    Animal ID from the session folder name
  extractEpochPower.m              Power extraction per epoch
  storeRunPower.m / storeNonRunPower.m
  buildRow.m                       Builds table rows from epoch power
  robust_zscore_matrix.m           Robust z-score (median/MAD)
  removeLineNoise.m                Spectral-interpolation line-noise removal (50/100/150 Hz default)
  rippleMasterDetector.m           Ripple + sharp-wave detection wrapper
  findRipples_CNN_modpeak.m        CNN ripple detection (modified peak)
  computeRippleStats_new.m         Per-ripple descriptive statistics
  computeLocalizedMedian.m
  buzcode_eacheventCSD.m           Event-triggered CSD (per event)
  ps1to120_Paria.m                 Power spectrum 1–120 Hz per layer (used by the pipeline)
  ps1to120_Alberto.m               Variant of the above (used by analyzepsd_Alberto.m)
  ps1to120_Rut.m                   Older version (not called by anything)
  saveFigSession.m / saveGroupFig.m / savePSFig.m   Figure helpers

optional_analysis/                 Stand-alone analyses, run after the main pipeline
  edit_hippocampalLayers_alberto.m           Edit layer channel assignments
  rippleanalysis/
    analyzecsdripples_Alberto.m              CSD of ripples + mean waveforms per recording
    analyzerippleslope_Alberto.m             Ripple slope check (pyramidal-channel position bias)
    ripplefinalscript.m                      Ripple UMAP pre-processing and analysis
  spectrumanalysis/
    analyzepsd_Alberto.m                     PSD by group and layer (dB and µV²/Hz)
    analyzecfcthetagamma_Alberto.m           Theta–gamma cross-frequency coupling (Modulation Index)
  spikeanalysis/
    analyzecellmetrics_Alberto.m             Spike-metrics pipeline (cell counts, INH/EXC ratios, LME)
    mandatory_functions/                     Helpers for the analyses above (stats, outliers, CSD)

old_unused_scripts_or_functions/   Legacy scripts, kept for reference only
```

---

## Requirements

### Software

- MATLAB (version used: **TODO – fill in**)
- Toolboxes (likely needed — confirm): Signal Processing, Statistics and Machine Learning (`fitlme`/`fitglme`), Parallel Computing (only if you re-enable `parfor`)
- Python environment with the CNN ripple detector (`cnn-env`) and the [cnn-matlab](https://github.com/PridaLab/cnn-matlab) repository (PridaLab)

### External code (not included in this repo)

Add all of these to the MATLAB path before running:

- **[HippoCookBook](https://github.com/valegarman/HippoCookBook)** (Manu's code) — produces the input files listed below. Add it, with its subfolders, to the MATLAB path. It also provides `getLFP`, `bz_Filter`, `bz_Diff`, `computeCSD`, `run_umap`, `SleepScoreMaster` and other helpers used here.
- **[CellExplorer](https://github.com/petersenpeter/CellExplorer)** (a HippoCookBook dependency; add it with its subfolders) — provides `InIntervals`, `SubtractIntervals`, `ConsolidateIntervals`, `ProcessCellMetrics`, `loadSession` and `CCG`. Note that these are **not** in HippoCookBook itself.
- **[cnn-matlab](https://github.com/PridaLab/cnn-matlab)** (PridaLab) — provides `detect_ripples_cnn` and `get_intervals`, used by `findRipples_CNN_modpeak.m` and `rippleMasterDetector.m`.

### Paths to edit

In **P3.0** of `LFPpipeline.m`, set your own paths (to the local `cnn-matlab` copy and the Python executable of `cnn-env`):

```matlab
CNN_MODEL_PATH = 'C:\Users\usuario\cnn-matlab-master\cnn';
CNN_PYTHON_EXE = 'C:\Users\usuario\anaconda3\envs\cnn-env\python.exe';
```

---

## Expected input data

One folder per experimental group; inside each, **one subfolder per recording/session**:

```
GroupFolder/
├── 2026-06-09_16-47-11-F016/      ← session folder (name ends in -<AnimalID>)
│   ├── <rec>.ripples.events.mat
│   ├── <rec>.thetaEpochs6to12.states.mat
│   ├── <rec>_hippocampalLayers_channelinfo.mat   (fields: oriens, pyramidal, radiatum, slm; NaN = unassigned)
│   └── ... (outputs from Manu's code)
└── ...
```

- The **animal ID** is the last chunk after the final `-` in the session folder name (e.g. `2026-06-09_16-47-11-F016` → `F016`). It is used for the hierarchical outlier removal (recordings nested in animals).
- Layers analysed: `or`, `pyr`, `rad`, `slm`. If a layer is missing, that session gets `NaN` for it and the rest are processed.

---

## How to run the main pipeline

1. Run Manu's code on all recordings.
2. Add this repo (with subfolders) and the external code to the MATLAB path.
3. Edit the CNN paths (see above).
4. Open `LFPpipeline.m` and run it **section by section, in order** (P1 → P3 → P1R → P2). Dialog boxes will ask for:
   - number of groups, a label for each group, and the folder of each group;
   - an output folder for each part (Part I, ripples, power spectrum).
5. Group labels in Part II must match the ones used by `ps1to120_Paria`.

> The commented line `%clearvars; close all; clc;` at the top is intentional so that variables are not wiped when you run sections one at a time.

---

## Outputs

| Part | Output |
|------|--------|
| P1R | `theta_allrec.xlsx`, `theta_run.xlsx`, `theta_nonrun.xlsx`, `locevent_table.xlsx`, plus `figures/` (per session and per group) |
| P3 | New `<rec>.ripples.events.mat` per recording, ripple tables (per-ripple measures with `subject` column), outlier-filtered tables, histograms per group |
| P2 | Excel tables of peak frequency, power-to-peak and mean power per band × condition × layer (one row per recording, `NaN` for missing data) and figures in `figures/ps/` |

Power-spectrum conditions: `allrec`, `loc` (locomotion), `nonloc`, `nonlocnonripple`, each in **dB** and **µV**.

---

## Key parameters

Defined at the top of each part of `LFPpipeline.m`.

**Part I – Theta/gamma (P1.0)**

| Parameter | Default |
|-----------|---------|
| `THETA_BAND` | 4–12 Hz |
| `LGAMMA_BAND` / `HGAMMA_BAND` | 20–60 / 60–100 Hz |
| `N_CYCLES` | 8 |
| `MIN_RUN_DUR` | 1 s (shorter running epochs are discarded) |
| `DETECT_OPTS` | theta bandpass 6–12 Hz, delta 1–3.5 Hz, `powerThreshold` = 1 |

**Part III – Ripples (P3.0)**

| Parameter | Default |
|-----------|---------|
| `RIPPLE_THRESHOLD` | [2 5] SD (min, max) |
| `RIPPLE_PASSBAND` | 80–240 Hz |
| `PEAKNORM_CUTOFF` | 50 (ripples with `PeakNormPow` above this are discarded) |
| `OUTLIER_THRESH_MAD` | 3.5 (robust z-score, MAD) |
| `OUTLIER_MIN_N_*` | 20 recordings / subjects / group |
| `OUTLIER_VARS` | `PeakAmp`, `PeakNormPow`, `Duration`, `PeakFq` (the first three log-transformed) |

**Part II – Power spectrum (P2.0)**

| Band | Range (Hz) |
|------|-----------|
| delta | 1–3.5 |
| theta | 6–12 |
| beta | 12–30 |
| low gamma | 30–60 |
| medium gamma | 60–90 |
| high gamma | 90–120 |

Each band/condition can be switched off with the `active` flag in `PS_BANDS` / `PS_CONDITIONS`.

---

## Optional analyses

Stand-alone scripts, run **after** the main pipeline. Each one asks for its own folders at the start.

| Script | Purpose |
|--------|---------|
| `edit_hippocampalLayers_alberto.m` | Loops over session folders and lets you edit the `oriens / pyramidal / radiatum / slm` channels in each `*_hippocampalLayers_channelinfo.mat` from the console (Enter = keep value). |
| `rippleanalysis/analyzecsdripples_Alberto.m` | CSD of ripples and mean ripple waveforms per recording. |
| `rippleanalysis/analyzerippleslope_Alberto.m` | Checks whether the ripple (pyramidal) channel sits at the same anatomical position across recordings, without re-running detection, and compares ripple measures between groups adjusting for that position (recordings nested in animals). Input: `all_ripples_original.xlsx`. |
| `rippleanalysis/ripplefinalscript.m` | UMAP pre-processing and analysis of ripples. |
| `spectrumanalysis/analyzepsd_Alberto.m` | PSD (0.1–120 Hz) by group and layer, in dB and µV²/Hz. |
| `spectrumanalysis/analyzecfcthetagamma_Alberto.m` | Theta-phase / gamma-amplitude coupling (Modulation Index, Tort et al. 2010), surrogate z-scores, preferred phase and comodulograms, per layer and run/non-run. |
| `spikeanalysis/analyzecellmetrics_Alberto.m` | Spike-metrics pipeline: ProcessCellMetrics recompute, cleaning, cell counts and INH/EXC ratios, batch LME. |

---

## Known issues / TODO

- [ ] Add the MATLAB version used and confirm the required toolboxes.
- [ ] `ripplefinalscript.m` calls names that do not exist in this repo: `buzcode_eventCSD` (the file is `buzcode_eacheventCSD.m`) and `analyze_ripple_measure_lme`, `analyze_ripple_measure_glme` and `analyze_ripple_measure_lme_figures` (the files are `analyze_ripple_measure_glme_generalized.m` and `analyze_ripple_measure_lme_figures_generalized.m`; there is no `analyze_ripple_measure_lme`). Update the calls or the file names.
- [ ] `rippleMasterDetector.m` has a `'filtercnn'` branch that calls `findRipples_wCNN`, which is not available anywhere. The pipeline uses the `findRipples_CNN_modpeak` branch, so it is not affected.
- [ ] Hard-coded CNN/Python paths in P3.0.
- [ ] Three versions of `ps1to120` (`_Paria`, `_Alberto`, `_Rut`): keep only the ones still in use (`_Rut` is not called by anything).
- [ ] `optional_analysis/spikeanalysis/mandatory_functions/` contains ripple/PSD statistics and scripts, not only spike helpers; consider renaming/splitting (e.g. `stats/`).
- [ ] `spike_analysis_pipeline.m` overlaps with `analyzecellmetrics_Alberto.m`; keep one.
- [ ] `old_unused_scripts_or_functions/` no longer holds anything the pipeline needs (`runThetaStage` now lives in `mandatory_functions/`); it can be removed since git keeps the history.

---

## Contact

L234 lab, CNC. For questions about the pipeline, open an issue in this repository or contact the lab members listed in the commit history.
