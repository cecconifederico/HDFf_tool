%% RUN for Cumulative Fréchet Distance 
%% (Saves HFD(K) for each channel up to kmax) 
clearvars
clc
close all

% --- PROTOCOL PARAMETERS ---
Fs = 200; 

kmax = 35;
kmin = 2;
K_range = kmin:kmax;

alpha = 2;      % Function parameter 
minPtsFit = 8;  % Minimum valid points to perform the linear fit 

datasetList = {'A1', 'M1', 'S1'};

HFDf_totale_1 = table();

%% Loop to load data: A1, M1, S1
for d = 1:length(datasetList)
    currentDS = datasetList{d};
    fprintf('\n--- Processing Dataset: %s ---\n', currentDS);

    msg = sprintf('Select the CSV file for %s', currentDS);
    [file, path] = uigetfile('*.csv', msg);
    if isequal(file,0)
        fprintf('Skipping dataset %s (selection canceled)\n', currentDS);
        continue;
    end

    data = readtable(fullfile(path, file));
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

        HFD_stack = [];   % Each row = HFD_vs_K of a valid segment
        nValidSeg = 0;

        for s = 1:nSeg
            segmento = segnale(startIdx(s):endIdx(s));
            segmento = segmento(:);
            segmento = segmento(isfinite(segmento));

            if numel(segmento) <= (2 * kmax)
                continue;
            end

            % --- Function call: obtain Lk ---
            [~, Lk, fitInfo] = fractaldim_frechet(segmento, kmax, alpha);

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

            HFD_stack = [HFD_stack; HFD_vs_K]; 
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

            res = table({currentDS}, {nomeAttuale}, {tipo}, {soggetto}, {elettrodo}, ...
                {HFD_median_curve}, {K_range}, ...
                'VariableNames', {'Area','Channel','ElectrodeType','SubjectID','Electrode','HFD_vs_K','K_range'});

            HFDf_totale_1 = [HFDf_totale_1; res]; 

            fprintf('Channel %s completed: HFD(K) curve saved (nValidSeg=%d)\n', ...
                nomeAttuale, nValidSeg);
        else
            fprintf('Channel %s: no valid HFD curve found\n', nomeAttuale);
        end
    end
end

save('HFDf_totale_curves.mat', 'HFDf_totale_1');