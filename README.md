# HDFf_tool – HFDncd

This branch contains MATLAB code implementing a Higuchi-based approach using the Normalized Compression Distance (NCD), referred to as HFDncd.

The method uses the Prediction by Partial Matching compressor (PPMd, variant H) to assess the similarity between an original signal and its low-pass-filtered and downsampled representations across different temporal scales.

Unlike HFDf, the method returns NCD values at each subsampling factor rather than a fractal dimension estimated by log-log regression.

The tool is designed for the analysis of stereo-electroencephalographic (sEEG) signals from different cortical areas and vigilance states. Both run scripts apply the same core computational procedure and the same within-state analysis to wakefulness and sleep recordings (N2, N3 and REM).

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
    ├── ncd_ppmd_bt.m
    ├── bw_lowpass_causal.m
    ├── run_HFDncd_awake.m
    └── run_HFDncd_sleep.m
```

## Main files

- `src/ncd_ppmd_bt.m`: main MATLAB function used to compute NCD across subsampling factors using the PPMd compressor.
- `src/bw_lowpass_causal.m`: MATLAB function implementing a causal Butterworth low-pass filter before subsampling.
- `src/run_HFDncd_awake.m`: script for wakefulness analysis, visualization and statistical comparisons.
- `src/run_HFDncd_sleep.m`: script applying the corresponding analysis separately to N2, N3 and REM recordings, including comparisons between sleep stages.
- `docs/input_format.md`: description of the expected input data format.
- `docs/pseudocode.md`: algorithmic specification and pseudocode.
- `docs/usage.md`: instructions for configuring and running the MATLAB scripts.
- `results/README.md`: description of the generated outputs.

## Data

The sEEG data used in this project are provided by the public database of the Montreal Neurological Institute.

The analysis focuses on three primary cortical areas:

- A1: Primary auditory cortex.
- M1: Primary motor cortex.
- S1: Primary somatosensory cortex.

Further details about the expected input format are available in [Input format](docs/input_format.md).

## Requirements

- MATLAB
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Python (a version supported by the MATLAB release in use)
- pyppmd Python library

Install pyppmd in the Python environment used by MATLAB:

```bash
python -m pip install pyppmd
```

## Usage

Add the `src` folder to the MATLAB path and configure the Python environment.

Run the sleep analysis first if the Awake–REM comparison is required by the wakefulness script.

The reference configuration uses `riportaBase = false` for all vigilance states. The downsampled signal is retained at its reduced sampling rate without restoring its original length.

Detailed instructions are available in [Usage](docs/usage.md).

## Output

The scripts save the channel-level NCD curves and analysis parameters as MATLAB `.mat` files:

- `NCDf_totale_curves.mat`: wakefulness results.
- `NCDppmdH_sleep_noBase.mat`: sleep-state results.

The scripts also generate plots and statistical tables in the MATLAB workspace.

Further information is available in [Results structure](results/README.md).

## Documentation

Detailed documentation is available in the `docs/` folder:

- [Input format](docs/input_format.md)
- [Algorithm and pseudocode](docs/pseudocode.md)
- [Usage](docs/usage.md)

The output data structures are described in [Results structure](results/README.md).

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

The implementation relies on the external Python library pyppmd, which is distributed under its own license.



The HFDncd implementation relies on the external Python library pyppmd, which is distributed under its own license.
