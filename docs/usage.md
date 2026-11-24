# Usage

This document explains how to run the HFDncd analysis in MATLAB using the Normalized Compression Distance (NCD) and the PPMd compressor.

## Requirements

The following software and libraries are required:

- MATLAB
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Python (compatible with the MATLAB version in use)
- pyppmd Python library

## Python Setup

Install the pyppmd library by executing the following command in the system terminal:

```bash
python -m pip install pyppmd
```

Open MATLAB and check the Python environment:

```matlab
pyenv
```

If necessary, configure MATLAB to use the Python executable corresponding to the environment in which pyppmd is installed:

```matlab
pyenv('Version', 'path_to_python_executable')
```

Verify that the pyppmd library can be imported correctly:

```matlab
py.importlib.import_module('pyppmd');
```

## MATLAB Setup

Open MATLAB and set the main repository folder `HDFf_tool` as the current working directory.

Add the `src` folder to the MATLAB path:

```matlab
addpath('src')
```

## Run the wakefulness analysis

To run the wakefulness analysis, execute:

```matlab
run_HFDncd_awake
```

When prompted, select the input files containing the wakefulness data for each cortical area.

The expected wakefulness input files are organized by cortical area:

```text
A1_rep.xlsx
M1_rep.xlsx
S1_rep.xlsx
```

The script supports both `.csv` and `.xlsx` input files.

Each file contains the recording channels associated with the corresponding cortical area.

## Run the sleep-state analysis

To run the sleep-state analysis, execute:

```matlab
run_HFDncd_sleep
```

When prompted, select the folder containing the sleep-state input files.

The script automatically searches for files according to cortical area and sleep stage.

The expected file naming convention is:

```text
A1_N2*.csv
A1_N3*.csv
A1_REM*.csv

M1_N2*.csv
M1_N3*.csv
M1_REM*.csv

S1_N2*.csv
S1_N3*.csv
S1_REM*.csv
```

Both `.csv` and `.xlsx` files are supported.

The wildcard character `*` allows additional characters in the file name after the sleep-stage identifier.

The script processes all matching files for each combination of cortical area and sleep stage.

The sleep stages considered are:

- N2: Non-REM sleep stage 2.
- N3: Non-REM sleep stage 3.
- REM: Rapid Eye Movement sleep.

The expected input format is described in:

```text
docs/input_format.md
```

## Parameters

The main parameters used for both wakefulness and sleep-state analyses are:

```matlab
Fs = 200;
ks = [2 5 10 20 33 67];
fc_Hz = [35 14 7 3.5 2 1];
n_bits = 3;
riportaBase = false;
minSegLen = 400;
k_stat = 2;
```

| Parameter | Description |
| --- | --- |
| `Fs` | Original sampling frequency (Hz). |
| `ks` | Subsampling factors. |
| `fc_Hz` | Butterworth low-pass cutoff frequency for each subsampling factor (Hz). |
| `n_bits` | Number of quantization bits. |
| `riportaBase` | Flag controlling restoration of the original temporal grid after subsampling. |
| `minSegLen` | Segment length threshold in samples. Only segments longer than this threshold are processed. |
| `k_stat` | Reference subsampling factor used for statistical analysis. |

The following computational parameters are defined in the core functions:

| Parameter | Value |
| --- | --- |
| Butterworth filter order | 7 |
| Filtering mode | Causal, single pass |
| Normalization | Z-score |
| Clipping threshold | +/- 4 standard deviations |
| Quantization levels | 8 |
| PPMd variant | H |
| Maximum context order | 6 |

These parameters can be modified directly in the MATLAB scripts and core functions if needed.

When `riportaBase = false`, the downsampled signal is retained at its reduced sampling rate without restoring the original temporal grid.

## Output

The scripts save the computed NCD curves and the corresponding analysis parameters as MATLAB `.mat` files in the current working directory.

### Wakefulness analysis

The output file is:

```text
NCDf_totale_curves.mat
```

The primary MATLAB table is:

```text
NCDf_totale_1
```

The table contains one row for each valid recording channel and seven columns describing the channel metadata and the corresponding NCD values.

### Sleep-state analysis

The output file is:

```text
NCDppmdH_sleep_noBase.mat
```

The primary MATLAB table is:

```text
NCDppmdH_sleep_noBase
```

The table contains one row for each valid recording channel processed from the input files and nine columns describing the channel metadata, sleep stage, source file, and corresponding NCD values.

For both wakefulness and sleep-state analyses, the NCD values are computed independently for each valid signal segment.

The median NCD value across valid segments is calculated for each subsampling factor and stored in the corresponding channel-level output.

Both output files also contain the following computational parameters:

```text
ks
fc_Hz
n_bits
riportaBase
Fs
```

Further details about the output data structures are available in:

```text
results/README.md
```

## Statistical Analysis and Visualization

The scripts generate graphical and statistical outputs for the investigation of neurodynamic differences across cortical areas and vigilance states.

### Wakefulness analysis

The wakefulness analysis includes:

- Channel-level NCD curves.
- Subject-specific NCD curves.
- Mean NCD curves and standard deviations across channels.
- Within-subject scatter plots comparing cortical areas.
- Paired Wilcoxon signed-rank tests.
- Effect size estimates.
- Median and interquartile range (IQR).
- Exploratory statistical comparisons across all subsampling factors.

### Sleep-state analysis

The sleep-state analysis includes:

- Subject-level NCD tables for N2, N3, and REM.
- Mean NCD curves and standard deviations across subjects for each cortical area and sleep stage.
- Within-subject scatter plots comparing cortical areas.
- Channel-level median and interquartile range (IQR).
- One-tailed paired Wilcoxon signed-rank tests between cortical areas.
- One-tailed paired Wilcoxon signed-rank tests between sleep stages.

The cortical area comparisons are performed separately for each sleep stage:

- M1 > S1
- M1 > A1
- S1 > A1

The sleep-stage comparisons are performed separately for each cortical area:

- REM > N2
- N2 > N3

The reference subsampling factor for statistical analysis is `k = 2`.

Subject-level NCD values are obtained by averaging the NCD values across recording channels belonging to the same subject, cortical area, and vigilance state.

The statistical tables are generated in the MATLAB workspace and are not automatically saved in the primary output files.

### Awake-REM comparison

If the Awake-REM comparison section is retained in the wakefulness script, the sleep-state analysis must be executed first.

The file `NCDppmdH_sleep_noBase.mat` must be available in the current working directory.

The comparison is performed within subjects and separately for each cortical area using a one-tailed paired Wilcoxon signed-rank test.

Both analyses must use the same computational parameters.

## Algorithm

The main NCD computation is implemented in:

```text
src/ncd_ppmd_bt.m
```

The causal Butterworth low-pass filtering is implemented in:

```text
src/bw_lowpass_causal.m
```

The method computes the Normalized Compression Distance between the original signal and its low-pass-filtered and downsampled representations using the PPMd compressor.

Unlike HFDf, HFDncd does not involve log-log regression or direct estimation of fractal dimension.

The output consists of NCD values computed at different subsampling factors.

A detailed description of the algorithm and its pseudocode is available in:

```text
docs/pseudocode.md
```
