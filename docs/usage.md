# Usage

This document explains how to run the HFDf analysis in MATLAB.

## Setup

Open MATLAB and set the main repository folder `HDFf_tool` as the current working directory.

Add the `src` folder to the MATLAB path:

```matlab
addpath('src')
```

## Run the wakefulness analysis

To run the wakefulness analysis, execute:

```matlab
run_HFDf
```

When prompted, select the input files containing the wakefulness data.

The expected wakefulness input files are organized by cortical area:

```text
A1_rep.xlsx
M1_rep.xlsx
S1_rep.xlsx
```

## Run the sleep-state analysis

To run the sleep-state analysis, execute:

```matlab
run_HFDf_sleep
```

When prompted, select the folder containing the CSV files organized by cortical area and sleep state.

The expected input format is described in:

```text
docs/input_format.md
```

## Parameters

The main parameters used in the scripts are:

```matlab
Fs = 200;
kmin = 2;
kmax = 35;
alpha = 2;
minPtsFit = 8;
```

These parameters can be modified directly in the scripts if needed.

## Output

The scripts save the output as MATLAB `.mat` files.

For the wakefulness analysis, the output file is:

```text
HFDf_totale.mat
```

For the sleep-state analysis, the output file is:

```text
HFDf_sleep.mat
```

Each output file contains a MATLAB table with one row for each valid channel and the corresponding HFDf(K) values.

## Algorithm

The HFDf computation is implemented in:

```text
src/fractaldim_frechet.m
```

A description of the algorithm is available in:

```text
docs/algorithm.md
```
