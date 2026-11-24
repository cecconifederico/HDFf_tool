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

Open MATLAB and check the configured Python environment:

```matlab
pyenv
```

If necessary, configure MATLAB to use the Python executable corresponding to the environment in which pyppmd is installed:

```matlab
pyenv('Version', 'path_to_python_executable')
```

Replace the placeholder with the actual path to the Python executable.

Configure the Python environment before using Python in MATLAB. Restart MATLAB if necessary.

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

## Run the Wakefulness Analysis

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

Both `.csv` and `.xlsx` files are supported.

Each file contains the recording channels associated with the corresponding cortical area.

## Run the Sleep-State Analysis

To run the sleep-state analysis, execute:

```matlab
run_HFDncd_sleep
```

When prompted, select the folder containing the sleep-state input files.

The script automatically searches the selected folder for files according to cortical area and sleep stage.

The expected file naming convention is:

```text
A1_N2.csv
M1_N2.csv
S1_N2.csv

A1_N3.csv
M1_N3.csv
S1_N3.csv

A1_REM.csv
M1_REM.csv
S1_REM.csv
```

Both `.csv` and `.xlsx` files are supported.

Additional characters after the sleep-stage identifier are accepted, for example:

```text
A1_N2_rep.csv
```

All matching files in the selected folder are processed. Avoid including duplicate copies of the same recording channels.

Missing area-stage combinations are skipped with a warning.

The same computational procedure is applied separately to each sleep stage:

- N2: Non-REM sleep stage 2.
- N3: Non-REM sleep stage 3.
- REM: Rapid Eye Movement sleep.

The expected input format for both analyses is described in [Input format](input_format.md).

## Execution Order

The wakefulness script includes an Awake–REM statistical comparison and loads the sleep-state output file:

```text
HFDncd_sleep.mat
```

Therefore, when performing the complete analysis, execute the sleep-state script before the wakefulness script:

```matlab
run_HFDncd_sleep
run_HFDncd_awake
```

Both scripts should be executed using the same MATLAB working directory.

For wakefulness-only analysis, the Awake–REM comparison section can be disabled.

## Parameters

Both run scripts use the same reference computational parameters:

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
| `k_stat` | Reference subsampling factor used for statistical analysis. Must be included in `ks`. |

The following computational parameters are defined in the core functions:

| Parameter | Value |
| --- | --- |
| Butterworth filter order | 7 |
| Filtering mode | Causal, single pass |
| Normalization | Independent z-score |
| Clipping threshold | ±4 standard deviations |
| Quantization levels | 8 |
| PPMd variant | H |
| Maximum context order | 6 |

When `riportaBase = false`, the downsampled signal is retained at its reduced sampling rate without restoring the original temporal grid.

These parameters can be modified directly in the MATLAB scripts and core functions if needed.

The same computational configuration must be used for all vigilance states when comparing their NCD values.

If the parameters are changed, regenerate the corresponding output files before performing between-state comparisons.

## Output

The scripts save the computed NCD curves and the corresponding analysis parameters as MATLAB `.mat` files in the current working directory.

### Wakefulness Analysis

The output file is:

```text
HFDncd_awake.mat
```

The primary MATLAB table is:

```text
HFDncd_awake
```

The table contains one row for each retained recording channel and seven columns describing the channel metadata and the corresponding NCD values.

### Sleep-State Analysis

The output file is:

```text
HFDncd_sleep.mat
```

The primary MATLAB table is:

```text
HFDncd_sleep
```

The table contains one row for each retained recording channel and nine columns describing the channel metadata, sleep stage, source file and corresponding NCD values.

For both analyses, NCD values are computed independently for each valid signal segment.

The median NCD value across valid segments is calculated for each subsampling factor and stored in the corresponding channel-level output.

Both output files also contain the following computational parameters:

```text
ks
fc_Hz
n_bits
riportaBase
Fs
```

The output files can be loaded into MATLAB using:

```matlab
load('HFDncd_awake.mat');
load('HFDncd_sleep.mat');
```

Further information about the output data structures is available in [Results structure](../results/README.md).

## Statistical Analysis and Visualization

The same within-state statistical procedure is applied separately to Awake, N2, N3 and REM.

The analysis includes:

- Example channel-level and subject-specific NCD curves.
- Mean NCD curves and standard deviations across recording channels.
- Subject-level NCD values obtained by averaging channels within each subject and cortical area.
- Subject-level mean and standard deviation at the reference subsampling factor.
- Within-subject scatter plots comparing cortical areas.
- Channel-level median and interquartile range (IQR).
- Two-tailed paired Wilcoxon signed-rank tests using the approximate method.
- Effect size estimates.
- Exploratory statistical comparisons across all subsampling factors.

The cortical area comparisons are:

- S1 vs A1.
- M1 vs S1.
- M1 vs A1.

The effect size is calculated as:

$$
r=\frac{|Z|}{\sqrt{n}}
$$

The primary statistical comparisons are performed at `k_stat = 2`, when at least three paired observations are available.

Additional exploratory comparisons are performed for all subsampling factors without correction for multiple comparisons.

Subject-specific figures can optionally be exported as PNG files by setting:

```matlab
salvaFigure = true;
```

## Between-State Comparisons

Additional one-tailed paired Wilcoxon signed-rank tests are performed separately for each cortical area.

The directional comparisons are:

- Awake > REM.
- REM > N2.
- N2 > N3.

Subjects are matched across vigilance states using their anonymized subject identifiers.

The Awake–REM comparison requires the previously generated `HFDncd_sleep.mat` output file.

The statistical tables are generated in the MATLAB workspace and are not automatically saved in the primary output files.

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

A detailed description of the algorithm and its pseudocode is available in [Algorithmic specification](pseudocode.md).
