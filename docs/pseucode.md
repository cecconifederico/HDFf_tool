# Algorithmic Specification for HFDncd Core Functions

This reference block describes the mathematical procedure and
algorithmic steps of the Higuchi-based approach using the
Normalized Compression Distance (NCD), referred to as HFDncd.

The method evaluates the similarity between an original signal
and its low-pass-filtered and downsampled representations using
the Prediction by Partial Matching compressor (PPMd, variant H).

---

## 1. Main NCD Computation

```text
Procedure NCD_PPMD_BT(X, ks, n_bits, Fs, fc_Hz, riportaBase);
begin
    // 1. Inputs:
    //    X           (1D Signal Vector)
    //    ks          (Subsampling Factors)
    //    n_bits      (Number of Quantization Bits)
    //    Fs          (Original Sampling Frequency)
    //    fc_Hz       (Low-Pass Cutoff Frequencies)
    //    riportaBase (Temporal Grid Restoration Flag)

    // 2. Outputs:
    //    Lk          (NCD Array)
    //    info        (Computational Parameters)

    Initialize_Default_Parameters(
        ks := [2, 5, 10, 20, 33, 67],
        n_bits := 3,
        Fs := 200,
        fc_Hz := [35, 14, 7, 3.5, 2, 1],
        riportaBase := false
    );

    Initialize_Constants(
        filter_order := 7,
        variant := 'H',
        max_order := 6
    );

    If Length(fc_Hz) != Length(ks) Then
        RaiseError("Inconsistent scale and cutoff arrays");

    P := RemoveNonFiniteValues(X);
    N := Length(P);

    Initialize_Structure(info,
        N, ks, fc_Hz, Fs, n_bits,
        riportaBase, variant, max_order
    );

    Initialize_Array(Lk,
        size := Length(ks),
        value := NaN
    );

    If N < 20 Then
        Return (Lk, info);

    // Quantization and compression of the original signal

    Pq := LOCAL_QUANTIZE(P, n_bits);

    cP := LOCAL_PPMD_CSIZE(
        Pq, variant, max_order
    );

    info.cP := cP;

    // Main scale-loop for NCD computation

    For i := 1 To Length(ks) Do
    begin
        k := ks[i];
        fc := fc_Hz[i];

        If N < 2 * k Then
            Continue;

        // Anti-aliasing filtering

        Y := BW_LOWPASS_CAUSAL(
            P, fc, Fs, filter_order
        );

        // Signal subsampling

        q := Subsample(
            Y,
            from_index := 1,
            to_index := N,
            step := k
        );

        // Optional temporal grid restoration

        If riportaBase = true Then
        begin
            Qk := RepeatEachSample(q, k);
            Qk := Truncate(Qk, length := N);
        end
        Else
            Qk := q;

        // Quantization of the downsampled signal

        Qq := LOCAL_QUANTIZE(Qk, n_bits);

        // Compression of the individual signal
        // and concatenated sequences

        cQ := LOCAL_PPMD_CSIZE(
            Qq, variant, max_order
        );

        PQ := Concatenate(Pq, Qq);

        cPQ := LOCAL_PPMD_CSIZE(
            PQ, variant, max_order
        );

        // Normalized Compression Distance

        Lk[i] := (cPQ - Minimum(cP, cQ))
                 / Maximum(cP, cQ);

    end;

    Return (Lk, info);

end;
```

---

## 2. Causal Butterworth Low-Pass Filtering

```text
Procedure BW_LOWPASS_CAUSAL(X, fc, Fs, order);
begin
    // 1. Inputs:
    //    X     (1D Signal Vector)
    //    fc    (Cutoff Frequency)
    //    Fs    (Sampling Frequency)
    //    order (Butterworth Filter Order)

    // 2. Output:
    //    Y     (Filtered Signal)

    Initialize_Default_Parameters(order := 7);

    X := ConvertToDoubleColumnVector(X);

    // Cutoff frequency normalization

    Wn := fc / (Fs / 2);

    If (Wn <= 0) or (Wn >= 1) Then
        RaiseError("Invalid normalized cutoff frequency");

    // Butterworth filter design

    (z, p, gain) := Butterworth(
        order, Wn, type := "low"
    );

    // Conversion to second-order sections

    (SOS, g) := ZeroPoleToSOS(z, p, gain);

    Y := g * X;

    // Single-pass causal filtering

    For s := 1 To NumberOfRows(SOS) Do
    begin
        b := SOS[s, 1:3];
        a := SOS[s, 4:6];

        // Steady-state initial conditions

        y0 := Y[1];

        yss := y0 * Sum(b) / Sum(a);

        zi := ComputeFilterInitialConditions(
            b, a,
            past_outputs := [yss, yss],
            past_inputs := [y0, y0]
        );

        Y := ApplyCausalFilter(b, a, Y, zi);

    end;

    Return Y;

end;
```

---

## 3. Signal Normalization and Quantization

```text
Procedure LOCAL_QUANTIZE(X, n_bits);
begin
    // 1. Inputs:
    //    X      (1D Signal Vector)
    //    n_bits (Number of Quantization Bits)

    // 2. Output:
    //    B      (Quantized uint8 Sequence)

    Initialize_Constants(n_std := 4);

    X := ConvertToDoubleColumnVector(X);

    mu := Mean(X);
    sd := StandardDeviation(X);

    If sd = 0 Then
        sd := 1;

    // Z-score normalization

    Z := (X - mu) / sd;

    // Clipping at +/- 4 standard deviations

    Z := Maximum(Minimum(Z, n_std), -n_std);

    // Uniform quantization

    levels := 2 ^ n_bits;

    B := Round(
        ((Z + n_std) / (2 * n_std))
        * (levels - 1)
    );

    // Conversion to unsigned 8-bit integers

    B := ConvertToUInt8(B);

    Return B;

end;
```

---

## 4. PPMd Compression and Compressed Size Estimation

```text
Procedure LOCAL_PPMD_CSIZE(B, variant, max_order);
begin
    // 1. Inputs:
    //    B         (Quantized uint8 Sequence)
    //    variant   (PPMd Compressor Variant)
    //    max_order (Maximum Context Order)

    // 2. Output:
    //    C         (Compressed Size in Bytes)

    // Conversion from MATLAB to Python bytes

    PythonBytes := ConvertToPythonBytes(B);

    // Lossless compression using PPMd

    CompressedData := PPMd_Compress(
        PythonBytes,
        variant := variant,
        max_order := max_order
    );

    // Compressed sequence length

    C := Length(CompressedData);

    Return C;

end;
```

---

## 5. Mathematical Definition of NCD

For each subsampling factor k, the Normalized Compression
Distance between the original signal P and its processed
representation Qk is defined as:

$$
NCD(P,Q_k) =
\frac{C(PQ_k)-\min(C(P),C(Q_k))}
{\max(C(P),C(Q_k))}
$$

where:

- `C(P)` is the compressed size of the quantized original signal.
- `C(Qk)` is the compressed size of the quantized downsampled signal.
- `C(PQk)` is the compressed size of the concatenation of the two quantized signals.

All compressed sizes are expressed in bytes and obtained using
the PPMd compressor (variant H, maximum context order 6).

---

## 6. Reference Computational Parameters

| Parameter | Value |
| --- | --- |
| Original sampling frequency | 200 Hz |
| Subsampling factors | [2, 5, 10, 20, 33, 67] |
| Low-pass cutoff frequencies (Hz) | [35, 14, 7, 3.5, 2, 1] |
| Filter type | Causal Butterworth low-pass |
| Filter order | 7 |
| Normalization | Z-score |
| Clipping threshold | +/- 4 standard deviations |
| Quantization bits | 3 |
| Quantization levels | 8 |
| Compressor | PPMd-H |
| Maximum context order | 6 |
| Temporal grid restoration | false |
| Reference subsampling factor | 2 |

---

## 7. Methodological Notes

The downsampled signal is retained at its reduced sampling
rate without restoring the original temporal grid
(`riportaBase = false`).

The original and downsampled signals are independently
normalized and quantized before compression.

The same computational procedure is applied to recordings
from different cortical areas (A1, M1, and S1) and vigilance
states (Awake, N2, N3, and REM).

For each recording channel, NCD values are computed
independently for each valid signal segment. The median NCD
value across valid segments is then calculated for each
subsampling factor.

Unlike the HFDf approach, HFDncd does not involve log-log
linear regression or a direct estimation of fractal dimension.

The output is a vector of NCD values representing the
similarity between the original signal and its processed
representations across different temporal scales.
    end;
    
    d := CA[nP, nQ];
    Return d;
end;
