% RUN for Frechet-based HFDf analysis across sleep states %
clearvars
clc
close all

% --- PROTOCOL PARAMETERS ---
Fs = 200; %#ok<NASGU>

kmax = 35;
kmin = 2;
K_range = kmin:kmax;

alpha = 2;      % Function parameter
minPtsFit = 8;  % Minimum valid points to perform the linear fit

regionList = {'A1', 'M1', 'S1'};
stateList  = {'N2', 'N3', 'REM'};

HFDf_sleep = table();

% --- SELECT FOLDER CONTAINING CSV FILES ---
baseFolder = uigetdir(pwd, 'Select the folder containing the CSV files');

if isequal(baseFolder, 0)
    error('No folder selected.');
end

fprintf('\nSelected folder: %s\n', baseFolder);

% Loop to load data: A1, M1, S1 and sleep states
for d = 1:length(regionList)

    currentDS = regionList{d};

    for st = 1:length(stateList)

        currentState = stateList{st};

        fprintf('\n--- Processing Dataset: %s - %s ---\n', currentDS, currentState);

        % Expected file name, for example:
        % A1_N2_Representative.csv
        pattern = sprintf('%s_%s*.csv', currentDS, currentState);
        fileInfo = dir(fullfile(baseFolder, pattern));

        if isempty(fileInfo)
            fprintf('Skipping dataset %s - %s: file not found\n', currentDS, currentState);
            continue;
        end

        filePath = fullfile(baseFolder, fileInfo(1).name);
        fprintf('Selected file: %s\n', fileInfo(1).name);

        data = readtable(filePath);
        nomiCanali = data.Properties.VariableNames;

        for c = 1:length(nomiCanali)

            nomeAttuale = nomiCanali{c};
            segnale = data.(nomeAttuale);
            x = segnale(:);

            % --- Robust segmentation: gap = padding value (absolute minimum) ---
            gapVal = min(abs(x(isfinite(x))));

            if isempty(gapVal) || ~isfinite(gapVal)
                fprintf('Channel %s: invalid gapVal -> skipping\n', nomeAttuale);
                continue;
            end

            tol = 1.05;
            isGapRaw = abs(x) <= tol * gapVal;

            gapMin = 200;
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
            fprintf('Channel %s -> Segments found: %d\n', nomeAttuale, nSeg);

            HFD_stack = [];
            nValidSeg = 0;

            for s = 1:nSeg

                segmento = segnale(startIdx(s):endIdx(s));
                segmento = segmento(:);
                segmento = segmento(isfinite(segmento));

                if numel(segmento) <= (2 * kmax)
                    continue;
                end

                % --- Function call: obtain Lk ---
                [~, Lk, ~] = fractaldim_frechet(segmento, kmax, alpha);

                if isempty(Lk) || all(~isfinite(Lk))
                    continue;
                end

                % --- Build HFD_vs_K from Lk vector ---
                HFD_vs_K = NaN(1, numel(K_range));

                for Kfit = kmin:kmax

                    kk = (kmin:Kfit)';
                    yy = Lk(kk)';

                    valid = isfinite(yy) & (yy > 0);

                    if nnz(valid) < minPtsFit
                        continue;
                    end

                    xfit = log(1 ./ kk(valid));
                    yfit = log(yy(valid));

                    coeff = polyfit(xfit, yfit, 1);
                    HFD_vs_K(Kfit - kmin + 1) = coeff(1);
                end

                if all(~isfinite(HFD_vs_K))
                    continue;
                end

                HFD_stack = [HFD_stack; HFD_vs_K]; %#ok<AGROW>
                nValidSeg = nValidSeg + 1;
            end

            if ~isempty(HFD_stack)

                HFD_median_curve = median(HFD_stack, 1, 'omitnan');

                % Metadata extraction (GD/NG, Subject, Electrode)
                if numel(nomeAttuale) >= 7
                    tipo = nomeAttuale(1:2);
                    soggetto = nomeAttuale(3:5);
                    elettrodo = nomeAttuale(6:end-1);
                else
                    tipo = "";
                    soggetto = "";
                    elettrodo = "";
                end

                res = table( ...
                    {currentDS}, ...
                    {currentState}, ...
                    {nomeAttuale}, ...
                    {tipo}, ...
                    {soggetto}, ...
                    {elettrodo}, ...
                    {HFD_median_curve}, ...
                    {K_range}, ...
                    'VariableNames', ...
                    {'Area','State','Channel','ElectrodeType','SubjectID','Electrode','HFD_vs_K','K_range'});

                HFDf_sleep = [HFDf_sleep; res]; %#ok<AGROW>

                fprintf('Channel %s completed: HFD(K) curve saved (nValidSeg=%d)\n', ...
                    nomeAttuale, nValidSeg);
            else
                fprintf('Channel %s: no valid HFD curve found\n', nomeAttuale);
            end
        end
    end
end

save('HFDf_sleep.mat', 'HFDf_sleep');
