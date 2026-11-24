# HDFncd_tool


This branch contains MATLAB code implementing a Higuchi-based approach using the Normalized Compression Distance (NCD), referred to as HFDncd.

The method uses the Prediction by Partial Matching (PPMd) compressor, variant H, to assess the similarity between the original signal and its low-pass-filtered and downsampled versions, investigating the preservation of temporal organization across different scales.

The tool is designed for the analysis of sEEG signals from different cortical areas. The main function implements the NCD computation, while the run script applies the method to wakefulness data.

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
    ├── ncd_ppmd.m
    ├── bw_lowpass_causal.m
    └── run_HFDncd_awake.m
```

## Main files

- `src/ncd_ppmd_bt.m`: main MATLAB function used to compute NCD across different subsampling factors using the PPMd compressor.
- `src/bw_lowpass_causal.m`: MATLAB function implementing a causal Butterworth low-pass filter for anti-aliasing before subsampling.
- `src/run_HFDncd_awake.m`: script used to run the HFDncd analysis on wakefulness data, generate plots and perform statistical comparisons between cortical areas.
- `docs/input_format.md`: description of the expected input data format.
- `docs/pseudocode.md`: pseudocode and description of the HFDncd computation.
- `docs/usage.md`: instructions for running the MATLAB code.
- `results/README.md`: information about the output files.

## Data

The sEEG data used in this project are provided by the public database of the Montreal Neurological Institute.

The analysis focuses on three primary cortical areas:

- A1: Primary auditory cortex.
- M1: Primary motor cortex.
- S1: Primary somatosensory cortex.

Further details about the expected input format are available in:

```text
docs/input_format.md
```

## Requirements

- MATLAB
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Python (compatible with the MATLAB version in use)
- pyppmd Python library

The pyppmd library can be installed using:

```bash
pip install pyppmd
```

MATLAB must be configured to use the Python environment in which pyppmd is installed.

## Usage

The MATLAB code can be run after adding the `src` folder to the MATLAB path and configuring the Python environment.

Detailed instructions are available in:

```text
docs/usage.md
```

## Output

The run script saves the computed NCD curves and the corresponding analysis parameters in a MATLAB `.mat` file:

```text
NCDf_totale_curves.mat
```

The script also generates plots and performs statistical comparisons between cortical areas.

Further information about the generated outputs is available in:

```text
results/README.md
```

## Documentation

Detailed documentation is available in the `docs/` folder:

- `input_format.md`: describes the expected input data format.
- `pseudocode.md`: describes the HFDncd algorithm and its computational steps.
- `usage.md`: provides instructions for running the code.

## License

This project is licensed under the MIT License.

The HFDncd implementation relies on the external Python library pyppmd, which is distributed under its own license.
