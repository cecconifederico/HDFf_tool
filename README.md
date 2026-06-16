# HDFf_tool

This repository contains MATLAB code for the computation of a Higuchi-like Fractal Dimension modified with the discrete Frechet distance, referred to as HFDf.

The tool is designed for the analysis of sEEG signals from different cortical areas and physiological states. The main function implements the HFDf computation, while the run scripts apply the method to wakefulness and sleep-state datasets.

## Repository structure

```text
HDFf_tool/
├── README.md
├── LICENSE
├── docs/
│   ├── input_format.md
│   ├── pseudocode.md
│   └── usage.md
├── results/
│   └── README.md
└── src/
    ├── fractaldim_frechet.m
    ├── run_HFDf_awake.m
    └── run_HFDf_sleep.m
```

## Main files

- `src/fractaldim_frechet.m`: main MATLAB function used to compute HFDf.
- `src/run_HFDf_awake.m`: script used to run the HFDf analysis on wakefulness data.
- `src/run_HFDf_sleep.m`: script used to run the HFDf analysis on sleep-state data.
- `docs/input_format.md`: description of the expected input data format.
- `docs/pseudocode.md`: pseudocode and description of the HFDf computation.
- `docs/usage.md`: instructions for running the MATLAB scripts.
- `results/README.md`: information about the output files.

## Data

The sEEG data used in this project are provided by the public database of the Montreal Neurological Institute.

Further details about the expected input format are available in:

```text
docs/input_format.md
```

## Usage

The MATLAB scripts can be run after adding the `src` folder to the MATLAB path.

Detailed instructions are available in:

```text
docs/usage.md
```

## Output

The scripts save the results as MATLAB `.mat` files.

Further information about the generated outputs is available in:

```text
results/README.md
```

## Documentation

Detailed documentation is available in the `docs/` folder:

- `input_format.md`: describes the expected input data format.
- `pseudocode.md`: describes the HFDf algorithm and pseudocode.
- `usage.md`: explains how to run the code.

## License

This project is released under the MIT License.
