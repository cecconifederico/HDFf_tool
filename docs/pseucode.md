# Algorithmic Specification for HFDf Core Functions

This reference block defines the mathematical procedure and algorithmic steps for estimating the Higuchi Fractal Dimension via discrete Fréchet distance (HFDf)

---

```pascal
Procedure FRACTAL_DIM_FRECHET(X, kmax, alpha);
begin
    // 1. Inputs: X (1D Signal Vector), kmax (Max Scale), alpha (Scaling Exponent)
    // 2. Outputs: hfd (Fractal Dimension), Lk (Distance Array)
    
    Initialize_Constants(kmin := 2, minPtsFit := 8, R2min := 0.90);
    
    P := RemoveNonFiniteValues(X);
    N := Length(P);
    
    if N < 20 Then
    begin
        hfd := NaN;
        Lk := EmptyArray();
        Return (hfd, Lk);
    end;
    
    kmax := Minimum(kmax, Floor(N / 2));
    Initialize_Array(Lk, size := kmax, value := NaN);
    
    // Main scale-loop for L(k) computation
    For k := kmin To kmax Do
    begin
        Qk := Subsample(P, from_index := 1, to_index := N, step := k);
        
        If Length(Qk) < 2 Then
            Continue;
            
        dF := DISCRETE_FRECHET_DIST_1D(P, Qk);
        Lk[k] := dF / (k ^ alpha);
        
        If Lk[k] <= 0 Then
            Lk[k] := NaN;
    end;
    
    // Log-Log Linear Regression Fitting
    valid_indices := FindIndices(Lk = finite and Lk > 0);
    
    If Length(valid_indices) < minPtsFit Then
    begin
        hfd := NaN;
        Return (hfd, Lk);
    end;
    
    Initialize_Empty_Vectors(xfit, yfit);
    For each idx In valid_indices Do
    begin
        x_val := ln(1 / idx);
        y_val := ln(Lk[idx]);
        Append(xfit, x_val);
        Append(yfit, y_val);
    end;
    
    // Least-Squares Linear Fit: y = M * x + C
    (M, C) := Polyfit(xfit, yfit, degree := 1);
    hfd := M; 
    
    // Statistical Validation Rule
    R2 := ComputeCoefficientOfDetermination(xfit, yfit, M, C);
    If R2 < R2min Then
        hfd := NaN;
        
    Return (hfd, Lk);
end;


Procedure DISCRETE_FRECHET_DIST_1D(P, Q);
begin
    // 1. Inputs: P (Reference Curve), Q (Subsampled Curve)
    // 2. Output: d (Discrete Fréchet coupling distance)
    
    nP := Length(P);
    nQ := Length(Q);
    
    If (nP = 0) or (nQ = 0) Then
    begin
        d := 0;
        Return d;
    end;
    
    Initialize_Matrix(CA, rows := nP, cols := nQ, value := Infinity);
    CA[1, 1] := AbsoluteValue(P[1] - Q[1]);
    
    // Initialize first column boundary conditions
    For i := 2 To nP Do
        CA[i, 1] := Maximum(CA[i - 1, 1], AbsoluteValue(P[i] - Q[1]));
        
    // Initialize first row boundary conditions
    For j := 2 To nQ Do
        CA[1, j] := Maximum(CA[1, j - 1], AbsoluteValue(P[1] - Q[j]));
        
    // Dynamic Programming Matrix Traversal Space
    For i := 2 To nP Do
    begin
        For j := 2 To nQ Do
        begin
            dij := AbsoluteValue(P[i] - Q[j]);
            min_backtrack := Minimum(CA[i - 1, j], CA[i - 1, j - 1], CA[i, j - 1]);
            CA[i, j] := Maximum(min_backtrack, dij);
        end;
    end;
    
    d := CA[nP, nQ];
    Return d;
end;
