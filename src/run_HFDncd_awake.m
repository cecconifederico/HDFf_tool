%% RUN: NCD with PPMd-H on anti-aliased, decimated signals (awake)
% For each segment and each decimation factor k in ks: causal Butterworth
% low-pass filter (cutoff fc_Hz) -> downsampling by k -> NCD(P,Q) using PPMd-H.
% Requires: ncd_ppmd_bt.m and bw_lowpass_causal.m in the same folder,
%           a pyenv with Python 3.9-3.12 and pyppmd installed.
clearvars
clc
close all

%% --- PARAMETERS ---
Fs          = 200;                % Sampling frequency (Hz); EEG already low-passed at 80 Hz
ks          = [2 5 10 20 33 67];  % Decimation factors
fc_Hz       = [35 14 7 3.5 2 1];  % Butterworth cutoff for each k (Hz)
n_bits      = 3;                  % Quantization: 2^3 = 8 amplitude levels
reportBase  = false;              
minSegLen   = 400;                % Minimum segment length (samples)

K_range     = ks;
datasetList = {'A1', 'M1', 'S1'};
NCD_awake   = table();
fileUsed    = struct();

%% --- DATA LOADING AND NCD COMPUTATION ---
for d = 1:length(datasetList)
    currentDS = datasetList{d};
    fprintf('\n--- Processing Dataset: %s ---\n', currentDS);

    msg = sprintf('Select the file for %s', currentDS);
    [file, path] = uigetfile({'*.csv;*.xlsx', 'Data files (*.csv, *.xlsx)'}, msg);
    if isequal(file, 0)
        fprintf('Skipping dataset %s (selection canceled)\n', currentDS);
        continue;
    end
    fileUsed.(currentDS) = file;

    data = readtable(fullfile(path, file));
    channelNames = data.Properties.VariableNames;

    for c = 1:length(channelNames)
        currentName = channelNames{c};
        signal = data.(currentName);
        x = signal(:);

        % --- Segmentation: gaps are runs of the padding value (absolute minimum) ---
        gapVal = min(abs(x(isfinite(x))));
        if isempty(gapVal) || ~isfinite(gapVal)
            fprintf('Channel %s: invalid gapVal -> skipping\n', currentName);
            continue;
        end

        tol = 1.05;
        isGapRaw = abs(x) <= tol * gapVal;

        gapMin = 200;                       % Runs shorter than this are not gaps
        isGap = isGapRaw;

        d_gap = diff([0; isGapRaw; 0]);
        gs = find(d_gap == 1);
        ge = find(d_gap == -1) - 1;

        for ii = 1:numel(gs)
            if (ge(ii) - gs(ii) + 1) < gapMin
                isGap(gs(ii):ge(ii)) = false;
            end
        end

        isSignal = ~isGap & isfinite(x);
        d_idx = diff([0; isSignal; 0]);
        startIdx = find(d_idx == 1);
        endIdx   = find(d_idx == -1) - 1;

        nSeg = numel(startIdx);
        fprintf('Channel %s -> Segments found: %d\n', currentName, nSeg);

        NCD_stack = [];                     % Each row = NCD vs k of a valid segment
        nValidSeg = 0;

        for s = 1:nSeg
            segment = signal(startIdx(s):endIdx(s));
            segment = segment(:);
            segment = segment(isfinite(segment));

            if numel(segment) <= minSegLen
                continue;
            end

            % NCD for each k (causal filter + decimation + PPMd-H)
            [Lk, ~] = ncd_ppmd_bt(segment, ks, n_bits, Fs, fc_Hz, reportBase);

            if isempty(Lk) || all(~isfinite(Lk))
                continue;
            end

            NCD_stack = [NCD_stack; Lk]; %#ok<AGROW>
            nValidSeg = nValidSeg + 1;
        end

        if ~isempty(NCD_stack)
            % Median across segments -> one NCD(k) curve per channel
            NCD_median_curve = median(NCD_stack, 1, 'omitnan');

            % Metadata extraction (electrode type, subject, electrode)
            if numel(currentName) >= 7
                elType    = currentName(1:2);
                subject   = currentName(3:5);
                electrode = currentName(6:end-1);
            else
                elType = ''; subject = ''; electrode = '';
            end

            res = table({currentDS}, {currentName}, {elType}, {subject}, {electrode}, ...
                {NCD_median_curve}, {K_range}, ...
                'VariableNames', {'Area','Channel','ElectrodeType','SubjectID','Electrode','NCD_vs_K','K_range'});
            NCD_awake = [NCD_awake; res]; %#ok<AGROW>

            fprintf('Channel %s completed: NCD(k) curve saved (nValidSeg=%d)\n', ...
                currentName, nValidSeg);
        else
            fprintf('Channel %s: no valid NCD curve found\n', currentName);
        end
    end
end

% --- Sanity check: the same file must not be selected for two areas ---
fn = fieldnames(fileUsed);
for ii = 1:numel(fn)
    for jj = ii+1:numel(fn)
        if strcmp(fileUsed.(fn{ii}), fileUsed.(fn{jj}))
            warning('Same file selected for %s and %s (%s). Please check.', ...
                fn{ii}, fn{jj}, fileUsed.(fn{ii}));
        end
    end
end

if isempty(NCD_awake)
    error('No NCD curve computed: check input files and pyenv.');
end

save('NCD_awake_curves.mat', 'NCD_awake', 'ks', 'fc_Hz', 'n_bits', 'reportBase', 'Fs');
