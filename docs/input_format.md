# Input Format

This document describes the expected input format for the sEEG data used by the example analysis scripts.
The data used in this project are derived from the public MNI sEEG database from the Montreal Neurological Institute

## Data files

The analysis scripts expect input data organized by cortical area and physiological state.

Two types of input files are used:

1. wakefulness data, provided as Excel files (`.xlsx`);
2. sleep-state data, provided as CSV files (`.csv`).

## Wakefulness data

For the wakefulness analysis, the input files are organized by cortical area only.

Wakefulness data are provided as Excel files (`.xlsx`) using the following naming convention:

```text
<Area>_rep.xlsx
```

Examples:

```text
A1_rep.xlsx
M1_rep.xlsx
S1_rep.xlsx
```

Each file contains the sEEG signals recorded during wakefulness for one specific cortical area.

## Sleep-state data

For the sleep-state analysis, the input files are organized by cortical area and sleep state.

Sleep-state data are provided as CSV files (`.csv`) using the following naming convention:

```text
<Area>_<State>_Representative.csv
```

Examples:

```text
A1_N2_Representative.csv
A1_N3_Representative.csv
A1_REM_Representative.csv
M1_N2_Representative.csv
M1_N3_Representative.csv
M1_REM_Representative.csv
S1_N2_Representative.csv
S1_N3_Representative.csv
S1_REM_Representative.csv
```

Each file contains the sEEG signals recorded for one specific cortical area and one specific sleep state.

## Cortical areas

The expected cortical areas are:

* `A1`: primary auditory cortex;
* `M1`: primary motor cortex;
* `S1`: primary somatosensory cortex.

## Physiological states

The expected physiological states are:

* wakefulness;
* `N2`;
* `N3`;
* `REM`.

Wakefulness data are identified by the file name structure `<Area>_rep.xlsx`.

Sleep-state data are identified by the file name structure `<Area>_<State>_Representative.csv`.

## File structure

Each input file contains the sEEG signals related to a specific cortical area and physiological state.

Each column corresponds to one sEEG channel, and each row corresponds to one time sample.

In the current data organization, each subject is represented by one channel for each cortical area and physiological state.

Example:

| Channel_1 | Channel_2 | Channel_3 |
| --------: | --------: | --------: |
|     0.123 |     0.087 |     0.154 |
|     0.118 |     0.091 |     0.149 |
|     0.121 |     0.088 |     0.151 |

## Channel naming convention

The script extracts metadata from the channel name using the following structure:

```matlab
tipo = nomeAttuale(1:2);
soggetto = nomeAttuale(3:5);
elettrodo = nomeAttuale(6:end-1);
```

Therefore, the channel name should contain:

* electrode type in the first two characters;
* subject ID in characters 3 to 5;
* electrode label in the remaining characters.

Example:

```text
GD029Rd10N
```

For wakefulness data, channel names may contain a final `W`, for example:

```text
GD029Rd10W
```

## Gap handling

Some signals include additional padding or gap samples.

For continuous signals, the analysis is performed on the full valid signal segment.

For signals containing interruptions, the script detects and removes the gap regions before computing HFDf. The remaining valid signal portions are analyzed as separate segments.
