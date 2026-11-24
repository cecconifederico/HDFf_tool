%% ============================================================
% RUN SLEEP: NCD-PPMd
%
% Adattamento del codice AWAKE ai dati del SONNO
%
% Stati:
%   N2
%   N3
%   REM
%
% Aree:
%   A1
%   M1
%   S1
%
% Per ogni segmento e per ogni k:
%
%   Butterworth passa-basso
%          ↓
%   sottocampionamento
%          ↓
%   Qk NON riportato alla base di P
%          ↓
%   NCD(P,Qk) con PPMd
%
% riportaBase = FALSE
%
% Richiede:
%   ncd_ppmd_bt.m
%   bw_lowpass_causal.m
%   Python + pyppmd
%
%% ============================================================

clearvars
clc
close all


%% ============================================================
% PARAMETRI
%% ============================================================

Fs = 200;

ks = [2 5 10 20 33 67];

fc_Hz = [35 14 7 3.5 2 1];

n_bits = 3;

% ------------------------------------------------------------
% Qk rimane sottocampionato:
% NON viene riportato alla stessa lunghezza temporale di P
% ------------------------------------------------------------

riportaBase = false;

minSegLen = 400;

% k usato per scatterplot e statistiche
k_stat = 2;


K_range = ks;


idx_k = find(ks == k_stat,1);

if isempty(idx_k)

    error('k_stat deve essere uno dei valori presenti in ks.');

end


%% ============================================================
% AREE E STATI
%% ============================================================

areas = {'A1','M1','S1'};

states = {'N2','N3','REM'};


%% ============================================================
% SELEZIONE CARTELLA DEL SONNO
%% ============================================================

sleepFolder = uigetdir( ...
    pwd, ...
    'Seleziona la cartella contenente i file del sonno');


if isequal(sleepFolder,0)

    error('Nessuna cartella selezionata.');

end


fprintf('\nCartella selezionata:\n%s\n',sleepFolder);


%% ============================================================
% TABELLA RISULTATI
%% ============================================================

NCDppmdH_sleep_noBase = table();


%% ============================================================
% CARICAMENTO E CALCOLO
%% ============================================================

for st = 1:numel(states)

    currentState = states{st};


    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('                     %s\n',currentState);
    fprintf('====================================================\n');


    %% ========================================================
    % LOOP AREE
    %% ========================================================

    for a = 1:numel(areas)

        currentArea = areas{a};


        fprintf('\n--- %s - %s ---\n', ...
            currentState, ...
            currentArea);


        %% ====================================================
        % CERCA AUTOMATICAMENTE FILE
        %% ====================================================

        patternCSV = sprintf( ...
            '%s_%s*.csv', ...
            currentArea, ...
            currentState);


        patternXLSX = sprintf( ...
            '%s_%s*.xlsx', ...
            currentArea, ...
            currentState);


        filesCSV = dir( ...
            fullfile( ...
                sleepFolder, ...
                patternCSV));


        filesXLSX = dir( ...
            fullfile( ...
                sleepFolder, ...
                patternXLSX));


        files = [filesCSV; filesXLSX];


        %% ====================================================
        % NESSUN FILE
        %% ====================================================

        if isempty(files)

            warning( ...
                'Nessun file trovato per %s - %s', ...
                currentArea, ...
                currentState);

            continue

        end


        fprintf( ...
            'Trovati %d file per %s - %s\n', ...
            numel(files), ...
            currentArea, ...
            currentState);


        %% ====================================================
        % LOOP EVENTUALI FILE MULTIPLI
        %% ====================================================

        for f = 1:numel(files)

            fileName = files(f).name;


            fprintf('\nFile: %s\n',fileName);


            filePath = fullfile( ...
                sleepFolder, ...
                fileName);


            %% =================================================
            % LETTURA
            %% =================================================

            data = readtable(filePath);


            nomiCanali = ...
                data.Properties.VariableNames;


            %% =================================================
            % LOOP CANALI
            %% =================================================

            for c = 1:numel(nomiCanali)

                nomeAttuale = ...
                    nomiCanali{c};


                segnale = ...
                    data.(nomeAttuale);


                %% ---------------------------------------------
                % CONTROLLO NUMERICO
                %% ---------------------------------------------

                if ~isnumeric(segnale)

                    fprintf( ...
                        'Canale %s non numerico -> skip\n', ...
                        nomeAttuale);

                    continue

                end


                x = segnale(:);


                %% =================================================
                % CAMPIONI FINITI
                %% =================================================

                finiteValues = ...
                    x(isfinite(x));


                if isempty(finiteValues)

                    fprintf( ...
                        'Canale %s vuoto -> skip\n', ...
                        nomeAttuale);

                    continue

                end


                %% =================================================
                % IDENTIFICAZIONE GAP
                %% =================================================

                gapVal = ...
                    min(abs(finiteValues));


                if isempty(gapVal) || ~isfinite(gapVal)

                    fprintf( ...
                        'Canale %s: gapVal non valido -> skip\n', ...
                        nomeAttuale);

                    continue

                end


                tol = 1.05;


                isGapRaw = ...
                    abs(x) <= tol * gapVal;


                gapMin = 200;


                isGap = ...
                    isGapRaw;


                %% =================================================
                % GAP CONTIGUI
                %% =================================================

                d_gap = diff( ...
                    [0; isGapRaw; 0]);


                gs = find( ...
                    d_gap == 1);


                ge = find( ...
                    d_gap == -1) - 1;


                %% =================================================
                % GAP TROPPO CORTI:
                % NON vengono considerati veri gap
                %% =================================================

                for ii = 1:numel(gs)

                    if (ge(ii)-gs(ii)+1) < gapMin

                        isGap(gs(ii):ge(ii)) = false;

                    end

                end


                %% =================================================
                % SEGMENTI VALIDI
                %% =================================================

                isSignal = ...
                    ~isGap & isfinite(x);


                d_idx = diff( ...
                    [0; isSignal; 0]);


                startIdx = find( ...
                    d_idx == 1);


                endIdx = find( ...
                    d_idx == -1) - 1;


                nSeg = ...
                    numel(startIdx);


                fprintf( ...
                    'Canale %s -> Segmenti trovati: %d\n', ...
                    nomeAttuale, ...
                    nSeg);


                %% =================================================
                % MATRICE NCD SEGMENTI
                %% =================================================

                NCD_stack = [];

                nValidSeg = 0;


                %% =================================================
                % LOOP SEGMENTI
                %% =================================================

                for s = 1:nSeg

                    segmento = ...
                        segnale( ...
                            startIdx(s):endIdx(s));


                    segmento = segmento(:);


                    segmento = ...
                        segmento(isfinite(segmento));


                    %% ---------------------------------------------
                    % SEGMENTO TROPPO CORTO
                    %% ---------------------------------------------

                    if numel(segmento) <= minSegLen

                        continue

                    end


                    %% =================================================
                    % NCD PER TUTTI I k
                    %
                    % Butterworth + decimazione + PPMd
                    %
                    % riportaBase = false
                    %% =================================================

                    [Lk,~] = ...
                        ncd_ppmd_bt( ...
                            segmento, ...
                            ks, ...
                            n_bits, ...
                            Fs, ...
                            fc_Hz, ...
                            riportaBase);


                    %% ---------------------------------------------
                    % CONTROLLO
                    %% ---------------------------------------------

                    if isempty(Lk) || ...
                            all(~isfinite(Lk))

                        continue

                    end


                    %% ---------------------------------------------
                    % AGGIUNGE CURVA
                    %% ---------------------------------------------

                    NCD_stack = [ ...
                        NCD_stack; ...
                        Lk]; %#ok<AGROW>


                    nValidSeg = ...
                        nValidSeg + 1;

                end


                %% =================================================
                % MEDIANA TRA SEGMENTI
                %% =================================================

                if isempty(NCD_stack)

                    fprintf( ...
                        'Canale %s: nessuna curva NCD valida\n', ...
                        nomeAttuale);

                    continue

                end


                NCD_median_curve = ...
                    median( ...
                        NCD_stack, ...
                        1, ...
                        'omitnan');


                %% =================================================
                % METADATA DAL NOME DEL CANALE
                %% =================================================

                if numel(nomeAttuale) >= 7

                    tipo = ...
                        nomeAttuale(1:2);


                    soggetto = ...
                        nomeAttuale(3:5);


                    elettrodo = ...
                        nomeAttuale(6:end-1);

                else

                    tipo = '';

                    soggetto = '';

                    elettrodo = '';

                end


                %% =================================================
                % CREA RIGA TABELLA
                %% =================================================

                res = table( ...
                    {currentArea}, ...
                    {currentState}, ...
                    {fileName}, ...
                    {nomeAttuale}, ...
                    {tipo}, ...
                    {soggetto}, ...
                    {elettrodo}, ...
                    {NCD_median_curve}, ...
                    {K_range}, ...
                    'VariableNames', ...
                    { ...
                    'Area', ...
                    'State', ...
                    'File', ...
                    'Channel', ...
                    'ElectrodeType', ...
                    'SubjectID', ...
                    'Electrode', ...
                    'NCD_vs_K', ...
                    'K_range'});


                %% =================================================
                % AGGIUNGE ALLA TABELLA
                %% =================================================

                NCDppmdH_sleep_noBase = [ ...
                    NCDppmdH_sleep_noBase; ...
                    res]; %#ok<AGROW>


                fprintf( ...
                    ['Canale %s completato: curva NCD(k) ' ...
                     'salvata (nSegValidi=%d)\n'], ...
                    nomeAttuale, ...
                    nValidSeg);

            end

        end

    end

end


%% ============================================================
% CONTROLLO FINALE
%% ============================================================

if isempty(NCDppmdH_sleep_noBase)

    error([ ...
        'Nessuna curva NCD calcolata. ' ...
        'Controlla file, nomi, pyenv e pyppmd.']);

end


%% ============================================================
% SALVATAGGIO
%% ============================================================

save( ...
    'NCDppmdH_sleep_noBase.mat', ...
    'NCDppmdH_sleep_noBase', ...
    'ks', ...
    'fc_Hz', ...
    'n_bits', ...
    'riportaBase', ...
    'Fs');


fprintf('\n');
fprintf('====================================================\n');
fprintf('CALCOLO COMPLETATO\n');
fprintf('Canali totali: %d\n', ...
    height(NCDppmdH_sleep_noBase));
fprintf('File: NCDppmdH_sleep_noBase.mat\n');
fprintf('riportaBase = false\n');
fprintf('====================================================\n');


%% ============================================================
% COLORI
%% ============================================================

color_M1 = [0.85 0.15 0.15];

color_S1 = [0.15 0.50 0.80];

color_A1 = [0.10 0.45 0.20];


colors = struct( ...
    'A1',color_A1, ...
    'M1',color_M1, ...
    'S1',color_S1);


lineStyles = struct( ...
    'A1',':o', ...
    'M1','-s', ...
    'S1','--d');


%% ============================================================
% TABELLE PER SOGGETTO
% AL k DI RIFERIMENTO
%% ============================================================

subjectTables = struct();


for st = 1:numel(states)

    currentState = states{st};


    subjectTables.(currentState) = ...
        build_subject_table_sleep( ...
            NCDppmdH_sleep_noBase, ...
            currentState, ...
            idx_k);


    fprintf( ...
        '\n--- TABELLA SOGGETTI %s - k=%d ---\n', ...
        currentState, ...
        k_stat);


    disp( ...
        subjectTables.(currentState));

end


%% ============================================================
% CURVE MEDIE NCD(k)
%
% N2
% N3
% REM
%
% UNA SOPRA L'ALTRA
%% ============================================================

figure( ...
    'Color','w', ...
    'Position',[200 40 850 1050]);


tiledlayout( ...
    3, ...
    1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


for st = 1:numel(states)

    currentState = states{st};


    ax = nexttile;


    hold(ax,'on');

    grid(ax,'on');

    box(ax,'on');


    handles = [];

    labels = {};


    %% ========================================================
    % LOOP AREE
    %% ========================================================

    for a = 1:numel(areas)

        currentArea = areas{a};


        idx = ...
            strcmp( ...
                NCDppmdH_sleep_noBase.State, ...
                currentState) & ...
            strcmp( ...
                NCDppmdH_sleep_noBase.Area, ...
                currentArea);


        Ta = ...
            NCDppmdH_sleep_noBase(idx,:);


        if isempty(Ta)

            continue

        end


        %% ====================================================
        % IMPORTANTE:
        %
        % Qui calcoliamo prima una curva per soggetto,
        % facendo la MEDIA tra i canali dello stesso soggetto.
        %
        % Poi media +/- SD tra soggetti.
        %% ====================================================

        subjects = unique( ...
            Ta.SubjectID, ...
            'stable');


        subjectCurves = ...
            NaN( ...
                numel(subjects), ...
                numel(ks));


        for s = 1:numel(subjects)

            sid = subjects{s};


            idxSub = ...
                strcmp( ...
                    Ta.SubjectID, ...
                    sid);


            M = ...
                cell2mat( ...
                    Ta.NCD_vs_K(idxSub));


            subjectCurves(s,:) = ...
                mean( ...
                    M, ...
                    1, ...
                    'omitnan');

        end


        %% ====================================================
        % MEDIA +/- SD TRA SOGGETTI
        %% ====================================================

        m = ...
            mean( ...
                subjectCurves, ...
                1, ...
                'omitnan');


        sd = ...
            std( ...
                subjectCurves, ...
                0, ...
                1, ...
                'omitnan');


        %% ====================================================
        % PLOT
        %% ====================================================

        h = plot( ...
            ax, ...
            ks, ...
            m, ...
            lineStyles.(currentArea), ...
            'Color', ...
            colors.(currentArea), ...
            'LineWidth',1.8, ...
            'MarkerSize',6, ...
            'MarkerFaceColor','none');


        errorbar( ...
            ax, ...
            ks, ...
            m, ...
            sd, ...
            'Color', ...
            colors.(currentArea), ...
            'LineWidth',0.9, ...
            'CapSize',4, ...
            'LineStyle','none', ...
            'HandleVisibility','off');


        handles = [ ...
            handles, ...
            h]; %#ok<AGROW>


        labels = [ ...
            labels, ...
            currentArea]; %#ok<AGROW>

    end


    %% ========================================================
    % FORMATTAZIONE
    %% ========================================================

    title( ...
        ax, ...
        currentState, ...
        'FontWeight','bold', ...
        'FontSize',13);


    xlabel( ...
        ax, ...
        'k', ...
        'FontWeight','bold');


    ylabel( ...
        ax, ...
        'NCD', ...
        'FontWeight','bold');


    xticks(ax,ks);


    % piccolo spazio dopo k = 67
    xlim(ax,[0 72]);


    % limite superiore fisso
    ylimCurrent = ylim(ax);

    ylim(ax,[ylimCurrent(1) 1.06]);


    set( ...
        ax, ...
        'FontSize',11, ...
        'LineWidth',1.1, ...
        'GridAlpha',0.25);


    if ~isempty(handles)

        legend( ...
            ax, ...
            handles, ...
            labels, ...
            'Location','northeast');

    end

end


%% ============================================================
% SCATTERPLOT
%
% RIGHE:
%   N2
%   N3
%   REM
%
% COLONNE:
%   S1 vs A1
%   M1 vs S1
%   M1 vs A1
%
% k = k_stat
%% ============================================================

figure( ...
    'Name','Scatterplot NCD sleep', ...
    'Color','w', ...
    'Position',[80 50 1150 1000]);


confronti = { ...
    'A1','S1'; ...
    'S1','M1'; ...
    'A1','M1'};


limiti = [0.85 1.08];


colTri = { ...
    [0.55 0.82 0.92], ...
    [0.92 0.28 0.27], ...
    [0.92 0.28 0.27]};


for st = 1:numel(states)

    currentState = states{st};


    tab = ...
        subjectTables.(currentState);


    for p = 1:3

        areaX = ...
            confronti{p,1};


        areaY = ...
            confronti{p,2};


        X = ...
            tab.(areaX);


        Y = ...
            tab.(areaY);


        mask = ...
            ~isnan(X) & ...
            ~isnan(Y);


        X = X(mask);

        Y = Y(mask);


        ax = subplot( ...
            3, ...
            3, ...
            (st-1)*3+p);


        hold(ax,'on');


        if isempty(X)

            continue

        end


        %% ====================================================
        % TRIANGOLO Y > X
        %% ====================================================

        patch( ...
            ax, ...
            [limiti(1) limiti(2) limiti(1)], ...
            [limiti(1) limiti(2) limiti(2)], ...
            colTri{p}, ...
            'EdgeColor','none');


        %% ====================================================
        % DIAGONALE
        %% ====================================================

        plot( ...
            ax, ...
            limiti, ...
            limiti, ...
            '--', ...
            'Color',[0.7 0.7 0.7], ...
            'LineWidth',0.8);


        %% ====================================================
        % PUNTI
        %% ====================================================

        scatter( ...
            ax, ...
            X, ...
            Y, ...
            55, ...
            'filled', ...
            'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k');


        %% ====================================================
        % FORMAT
        %% ====================================================

        axis(ax,'square');

        box(ax,'on');


        xlim(ax,limiti);

        ylim(ax,limiti);


        xticks(ax,limiti);

        yticks(ax,limiti);


        xticklabels( ...
            ax, ...
            {'0.85','1.08'});


        yticklabels( ...
            ax, ...
            {'0.85','1.08'});


        xlabel( ...
            ax, ...
            areaX, ...
            'FontWeight','bold');


        ylabel( ...
            ax, ...
            areaY, ...
            'FontWeight','bold');


        if p == 2

            title( ...
                ax, ...
                currentState, ...
                'FontWeight','bold');

        end


        set( ...
            ax, ...
            'FontSize',10, ...
            'LineWidth',1, ...
            'Layer','top');

    end

end


%% ============================================================
% MEDIANA E IQR
%
% PER AREA E STATO
% AL k = k_stat
%
% SUI CANALI
%% ============================================================

riepilogo = table();


for st = 1:numel(states)

    currentState = states{st};


    for a = 1:numel(areas)

        currentArea = areas{a};


        idx = ...
            strcmp( ...
                NCDppmdH_sleep_noBase.State, ...
                currentState) & ...
            strcmp( ...
                NCDppmdH_sleep_noBase.Area, ...
                currentArea);


        Ta = ...
            NCDppmdH_sleep_noBase(idx,:);


        if isempty(Ta)

            continue

        end


        %% ====================================================
        % SOLO k = k_stat
        %% ====================================================

        vals = cellfun( ...
            @(x) x(idx_k), ...
            Ta.NCD_vs_K);


        vals = ...
            vals(isfinite(vals));


        if isempty(vals)

            continue

        end


        %% ====================================================
        % QUARTILI
        %% ====================================================

        q = ...
            prctile( ...
                vals, ...
                [25 50 75]);


        %% ====================================================
        % RIGA
        %% ====================================================

        riga = table( ...
            {currentState}, ...
            {currentArea}, ...
            numel(vals), ...
            q(2), ...
            q(1), ...
            q(3), ...
            iqr(vals), ...
            'VariableNames', ...
            { ...
            'State', ...
            'Area', ...
            'n_canali', ...
            'Mediana', ...
            'Q1', ...
            'Q3', ...
            'IQR'});


        riepilogo = [ ...
            riepilogo; ...
            riga]; %#ok<AGROW>

    end

end


fprintf('\n');
fprintf('====================================================\n');
fprintf('MEDIANA E IQR - k = %d\n',k_stat);
fprintf('====================================================\n');

disp(riepilogo);


%% ============================================================
% WILCOXON TRA AREE
%
% k = 2
%
% M1 > S1
% M1 > A1
% S1 > A1
%
% PAIRED
% ONE-TAILED
%% ============================================================

stats_sleep = table();


comparisons = { ...
    'M1','S1'; ...
    'M1','A1'; ...
    'S1','A1'};


for st = 1:numel(states)

    currentState = states{st};


    tab = ...
        subjectTables.(currentState);


    for c = 1:size(comparisons,1)

        areaHigh = ...
            comparisons{c,1};


        areaLow = ...
            comparisons{c,2};


        X = ...
            tab.(areaHigh);


        Y = ...
            tab.(areaLow);


        mask = ...
            ~isnan(X) & ...
            ~isnan(Y);


        X = X(mask);

        Y = Y(mask);


        n = ...
            numel(X);


        p = NaN;

        W = NaN;

        md = NaN;


        %% ====================================================
        % WILCOXON
        %% ====================================================

        if n >= 3

            [p,~,statW] = ...
                signrank( ...
                    X, ...
                    Y, ...
                    'tail','right');


            if isfield( ...
                    statW, ...
                    'signedrank')

                W = ...
                    statW.signedrank;

            end


            md = ...
                median(X-Y);

        end


        confronto = ...
            sprintf( ...
                '%s > %s', ...
                areaHigh, ...
                areaLow);


        riga = table( ...
            {currentState}, ...
            {confronto}, ...
            n, ...
            W, ...
            p, ...
            md, ...
            'VariableNames', ...
            { ...
            'State', ...
            'Comparison', ...
            'n', ...
            'SignedRank', ...
            'p_value', ...
            'median_difference'});


        stats_sleep = [ ...
            stats_sleep; ...
            riga]; %#ok<AGROW>

    end

end


fprintf('\n');
fprintf('====================================================\n');
fprintf('WILCOXON TRA AREE - k = %d\n',k_stat);
fprintf('====================================================\n');

disp(stats_sleep);


%% ============================================================
% WILCOXON TRA FASI DEL SONNO
%
% REM > N2
% N2  > N3
%
% PER:
% M1
% S1
% A1
%
% PAIRED
% ONE-TAILED
%% ============================================================

stageStats = table();


stageComparisons = { ...
    'REM','N2'; ...
    'N2','N3'};


for a = 1:numel(areas)

    currentArea = areas{a};


    for c = 1:size(stageComparisons,1)

        stateHigh = ...
            stageComparisons{c,1};


        stateLow = ...
            stageComparisons{c,2};


        tabHigh = ...
            subjectTables.(stateHigh);


        tabLow = ...
            subjectTables.(stateLow);


        %% ====================================================
        % ALLINEAMENTO SOGGETTI
        %% ====================================================

        [~,ia,ib] = ...
            intersect( ...
                tabHigh.SubjectID, ...
                tabLow.SubjectID, ...
                'stable');


        X = ...
            tabHigh.(currentArea)(ia);


        Y = ...
            tabLow.(currentArea)(ib);


        mask = ...
            ~isnan(X) & ...
            ~isnan(Y);


        X = X(mask);

        Y = Y(mask);


        n = ...
            numel(X);


        p = NaN;

        md = NaN;


        %% ====================================================
        % WILCOXON
        %% ====================================================

        if n >= 3

            p = ...
                signrank( ...
                    X, ...
                    Y, ...
                    'tail','right');


            md = ...
                median(X-Y);

        end


        confronto = ...
            sprintf( ...
                '%s > %s', ...
                stateHigh, ...
                stateLow);


        riga = table( ...
            {currentArea}, ...
            {confronto}, ...
            n, ...
            p, ...
            md, ...
            'VariableNames', ...
            { ...
            'Area', ...
            'Comparison', ...
            'n', ...
            'p_value', ...
            'median_difference'});


        stageStats = [ ...
            stageStats; ...
            riga]; %#ok<AGROW>

    end

end


fprintf('\n');
fprintf('====================================================\n');
fprintf('WILCOXON TRA FASI DEL SONNO - k = %d\n',k_stat);
fprintf('====================================================\n');

disp(stageStats);


%% ============================================================
% FUNZIONE LOCALE:
% TABELLA PER SOGGETTO
%% ============================================================

function tab = build_subject_table_sleep( ...
    T, ...
    stateName, ...
    idx_k)

% Costruisce:
%
% SubjectID | A1 | M1 | S1
%
% Per ogni soggetto:
% MEDIA tra i canali della stessa area
% al k selezionato.

    %% --------------------------------------------------------
    % SOLO STATO RICHIESTO
    %% --------------------------------------------------------

    idxState = ...
        strcmp( ...
            T.State, ...
            stateName);


    Ts = ...
        T(idxState,:);


    %% --------------------------------------------------------
    % SOGGETTI
    %% --------------------------------------------------------

    soggetti_unici = ...
        unique( ...
            Ts.SubjectID, ...
            'stable');


    tab = table( ...
        soggetti_unici, ...
        'VariableNames', ...
        {'SubjectID'});


    tab.A1 = ...
        NaN(height(tab),1);


    tab.M1 = ...
        NaN(height(tab),1);


    tab.S1 = ...
        NaN(height(tab),1);


    %% --------------------------------------------------------
    % LOOP SOGGETTI
    %% --------------------------------------------------------

    for i = 1:height(tab)

        sid = ...
            tab.SubjectID{i};


        %% ----------------------------------------------------
        % LOOP AREE
        %% ----------------------------------------------------

        for area = {'A1','M1','S1'}

            nm = ...
                area{1};


            idx = ...
                strcmp( ...
                    Ts.SubjectID, ...
                    sid) & ...
                strcmp( ...
                    Ts.Area, ...
                    nm);


            if any(idx)

                valori = cellfun( ...
                    @(x) x(idx_k), ...
                    Ts.NCD_vs_K(idx));


                tab.(nm)(i) = ...
                    mean( ...
                        valori, ...
                        'omitnan');

            end

        end

    end

end
