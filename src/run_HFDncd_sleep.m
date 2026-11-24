%% RUN SLEEP: NCD-PPMd-H (N2, N3, REM)
% Stessa pipeline del run_HFDncd_awake.m, applicata separatamente a ogni stato.
% Richiede ncd_ppmd_bt.m, bw_lowpass_causal.m e Python con pyppmd.
clearvars
clc
close all

%% --- PARAMETRI ---
Fs = 200;
ks = [2 5 10 20 33 67];
fc_Hz = [35 14 7 3.5 2 1];
n_bits = 3;
riportaBase = false;
minSegLen = 400;
k_stat = 2;
K_range = ks;
idx_k = find(ks == k_stat, 1);
if isempty(idx_k)
    error('k_stat deve essere uno dei valori di ks.');
end

areas = {'A1','M1','S1'};
states = {'N2','N3','REM'};
% Conteggi attesi nel dataset di riferimento (righe: stati; colonne: aree).
expectedCounts = [6 16 14; 6 16 14; 5 13 11];
salvaFigure = false;  % esportazione facoltativa PNG per soggetto

sleepFolder = uigetdir(pwd, 'Seleziona la cartella contenente i file del sonno');
if isequal(sleepFolder,0)
    error('Nessuna cartella selezionata.');
end

HFDncd_sleep = table('Size',[0 9], ...
    'VariableTypes',repmat({'cell'},1,9), ...
    'VariableNames',{'Area','State','File','Channel','ElectrodeType', ...
    'SubjectID','Electrode','NCD_vs_K','K_range'});

%% --- CARICAMENTO E CALCOLO (stessa segmentazione della veglia) ---
for st = 1:numel(states)
    currentState = states{st};
    for d = 1:numel(areas)
        currentDS = areas{d};
        fprintf('\n--- %s / %s ---\n', currentState, currentDS);
        files = [dir(fullfile(sleepFolder, [currentDS '_' currentState '*.csv'])); ...
                 dir(fullfile(sleepFolder, [currentDS '_' currentState '*.xlsx']))];
        if isempty(files)
            warning('Nessun file trovato per %s - %s.', currentDS, currentState);
            continue;
        end
        for f = 1:numel(files)
            fileName = files(f).name;
            data = readtable(fullfile(sleepFolder, fileName));
            nomiCanali = data.Properties.VariableNames;
            for c = 1:numel(nomiCanali)
                nomeAttuale = nomiCanali{c};
                segnale = data.(nomeAttuale);
                if ~isnumeric(segnale)
                    fprintf('Canale %s non numerico -> skip\n', nomeAttuale);
                    continue;
                end
                x = segnale(:);

                % --- segmentazione robusta: gap = valore di padding (minimo assoluto) ---
                gapVal = min(abs(x(isfinite(x))));
                if isempty(gapVal) || ~isfinite(gapVal)
                    fprintf('Canale %s: gapVal non valido -> skip\n', nomeAttuale);
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
                fprintf('Canale %s -> Segmenti trovati: %d\n', nomeAttuale, nSeg);

                NCD_stack = [];
                nValidSeg = 0;

                for s = 1:nSeg
                    segmento = segnale(startIdx(s):endIdx(s));
                    segmento = segmento(:);
                    segmento = segmento(isfinite(segmento));

                    if numel(segmento) <= minSegLen
                        continue;
                    end

                    % NCD per ogni k (filtro causale + decimazione + PPMd)
                    [Lk, ~] = ncd_ppmd_bt(segmento, ks, n_bits, Fs, fc_Hz, riportaBase);

                    if isempty(Lk) || all(~isfinite(Lk))
                        continue;
                    end

                    NCD_stack = [NCD_stack; Lk]; %#ok<AGROW>
                    nValidSeg = nValidSeg + 1;
                end

                if ~isempty(NCD_stack)
                    NCD_median_curve = median(NCD_stack, 1, 'omitnan');

                    if numel(nomeAttuale) >= 7
                        tipo = nomeAttuale(1:2);
                        soggetto = nomeAttuale(3:5);
                        elettrodo = nomeAttuale(6:end-1);
                    else
                        tipo = ''; soggetto = ''; elettrodo = '';
                    end
                    res = table({currentDS}, {currentState}, {fileName}, ...
                        {nomeAttuale}, {tipo}, {soggetto}, {elettrodo}, ...
                        {NCD_median_curve}, {K_range}, ...
                        'VariableNames', HFDncd_sleep.Properties.VariableNames);
                    HFDncd_sleep = [HFDncd_sleep; res]; %#ok<AGROW>
                    fprintf('Canale %s completato (nSegValidi=%d)\n', nomeAttuale, nValidSeg);
                else
                    fprintf('Canale %s: nessuna curva NCD valida\n', nomeAttuale);
                end
            end
        end
    end
end

if isempty(HFDncd_sleep)
    error('Nessuna curva NCD calcolata: controlla i file e pyenv.');
end
% Avviso per canali ripetuti nello stesso stato e area (es. CSV e XLSX duplicati).
[G,ga,gs,gc] = findgroups(HFDncd_sleep.Area, ...
    HFDncd_sleep.State, HFDncd_sleep.Channel);
counts = splitapply(@numel, HFDncd_sleep.Channel, G);
for i = find(counts > 1)'
    warning('Canale duplicato: %s / %s / %s. Controllare i file di origine.', ...
        ga{i}, gs{i}, gc{i});
end
save('NCDppmdH_sleep_noBase.mat', 'HFDncd_sleep', ...
    'ks', 'fc_Hz', 'n_bits', 'riportaBase', 'Fs');

%% --- ANALISI PER STATO: stesso codice della veglia ---
subjectTables = struct();
channelTables = struct();
stats_sleep = table();
tab_tutti_sleep = table();
riepilogo_sleep = table();
for st = 1:numel(states)
    statoLabel = states{st};
    fprintf('\n========== %s ==========\n', statoLabel);
    Tstate = HFDncd_sleep(strcmp(HFDncd_sleep.State, statoLabel), :);
    % Vista di lavoro con le stesse sette colonne e gli stessi nomi della veglia.
    % La tabella salvata mantiene State e File e lo schema originale del sonno.
    T = Tstate(:, {'Area','Channel','ElectrodeType','SubjectID','Electrode','NCD_vs_K','K_range'});
    T.Properties.VariableNames = {'Area','Canale','TipoElettrodo','SoggettoID','Elettrodo','NCD_vs_K','K_range'};
    channelTables.(statoLabel) = T;
    subjectTables.(statoLabel) = build_subject_table(T, idx_k);
    if isempty(T)
        warning('Nessun canale valido per %s: analisi saltata.', statoLabel);
        continue;
    end
    %% --- COLORI ---
    color_M1 = [0.85, 0.15, 0.15]; % Rosso opaco
    color_S1 = [0.15, 0.50, 0.80]; % Blu opaco
    color_A1 = [0.10, 0.45, 0.20]; % Verde scuro
    areas = {'A1', 'M1', 'S1'};

    %% --- PLOT ESEMPIO: primo canale ---
    figure; hold on; grid on; box on;
    plot(T.K_range{1}, T.NCD_vs_K{1}, 'o-', 'LineWidth', 2);
    xticks(ks);
    xlabel('k'); ylabel('NCD');
    title(sprintf('%s - Esempio curva NCD(k) - %s', statoLabel, T.Canale{1}));

    %% --- PLOT SOGGETTO SPECIFICO (mediana canali per area) ---
    subjectID = '029';
    if ~any(strcmp(T.SoggettoID, subjectID))
        subjectID = T.SoggettoID{1};
    end
    Tsub = T(strcmp(T.SoggettoID, subjectID), :);

    figure; hold on; grid on; box on
    plot_c = []; legend_labels = {};
    for i = 1:numel(areas)
        Ta = Tsub(strcmp(Tsub.Area, areas{i}), :);
        if isempty(Ta), continue; end
        Hmed = median(cell2mat(Ta.NCD_vs_K), 1, 'omitnan');
        switch areas{i}
            case 'A1', c = plot(ks, Hmed, '--o', 'Color', 'g', 'LineWidth', 2);
            case 'M1', c = plot(ks, Hmed, '-s',  'Color', 'r', 'LineWidth', 2);
            case 'S1', c = plot(ks, Hmed, ':d',  'Color', 'b', 'LineWidth', 2);
        end
        plot_c = [plot_c, c]; %#ok<AGROW>
        legend_labels = [legend_labels, areas{i}]; %#ok<AGROW>
    end
    legend(plot_c, legend_labels, 'Location', 'best');
    xticks(ks);
    xlabel('k', 'FontWeight', 'bold'); ylabel('NCD', 'FontWeight', 'bold');
    title(sprintf('%s - Sub %s', statoLabel, subjectID));

    %% --- TABELLA PER SOGGETTO AL k DI RIFERIMENTO ---
    tabella_soggetti = build_subject_table(T, idx_k);

    % controllo anti-duplicati (2): due aree identiche per tutti i soggetti?
    coppie_aree = {'A1','M1'; 'A1','S1'; 'M1','S1'};
    for ii = 1:size(coppie_aree, 1)
        v1 = tabella_soggetti.(coppie_aree{ii,1});
        v2 = tabella_soggetti.(coppie_aree{ii,2});
        mask = ~isnan(v1) & ~isnan(v2);
        if any(mask) && isequal(v1(mask), v2(mask))
            warning('I valori di %s e %s a k=%d sono IDENTICI per tutti i soggetti: probabile scambio di file.', ...
                coppie_aree{ii,1}, coppie_aree{ii,2}, k_stat);
        end
    end

    %% --- MEDIA +/- SD PER AREA (k = k_stat) ---
    areas_order = {'A1', 'S1', 'M1'};
    markers_fig = struct('A1', 'o', 'S1', 'd', 'M1', 's');
    colors_fig  = struct('A1', color_A1, 'S1', color_S1, 'M1', color_M1);

    figure('Color', 'w', 'Position', [100, 300, 420, 420]);
    hold on; box on; grid on;
    for i = 1:numel(areas_order)
        nm = areas_order{i};
        vals = tabella_soggetti.(nm);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        errorbar(i, mean(vals), std(vals), 'Color', colors_fig.(nm), 'LineWidth', 1.5, ...
            'CapSize', 10, 'LineStyle', 'none', 'HandleVisibility', 'off');
        plot(i, mean(vals), markers_fig.(nm), 'MarkerSize', 12, 'LineWidth', 1.5, ...
            'MarkerFaceColor', colors_fig.(nm), 'MarkerEdgeColor', colors_fig.(nm));
    end
    xlim([0.5, numel(areas_order) + 0.5]);
    xticks(1:numel(areas_order)); xticklabels(areas_order);
    ax = gca; ax.XAxis.FontWeight = 'bold';
    ylabel(sprintf('NCD (k=%d)', k_stat), 'FontWeight', 'bold', 'FontSize', 13);
    set(gca, 'FontSize', 11, 'LineWidth', 1.1);
    title(statoLabel);

    %% --- SCATTERPLOT CON DIAGONALE (k = k_stat) ---
    figure('Name', ['Scatterplots NCD: ' statoLabel], 'Color', 'w', 'Position', [100 300 1300 450]);

    confronti = {'A1', 'S1'; 'S1', 'M1'; 'A1', 'M1'};   % {asse X, asse Y}
    limiti    = [0.85 1.08];
    colTri    = {[0.55 0.82 0.92], [0.92 0.28 0.27], [0.92 0.28 0.27]};  % azzurro, rosso, rosso
    mkSize    = 70;     % dimensione pallini
    fsAsse    = 14;     % dimensione etichette assi
    fsTick    = 11;     % dimensione numeri sugli assi

    for p = 1:3
        areaX = confronti{p,1};
        areaY = confronti{p,2};
        X = tabella_soggetti.(areaX);
        Y = tabella_soggetti.(areaY);
        mask = ~isnan(X) & ~isnan(Y);
        X = X(mask); Y = Y(mask);

        subplot(1, 3, p); hold on;
        if isempty(X), continue; end

        % triangolo colorato sopra la diagonale (Y > X)
        patch([limiti(1) limiti(2) limiti(1)], ...
              [limiti(1) limiti(2) limiti(2)], ...
              colTri{p}, 'EdgeColor', 'none');

        % diagonale
        plot(limiti, limiti, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.75);

        % pallini neri
        scatter(X, Y, mkSize, 'filled', ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

        axis square; box on;
        xlim(limiti); ylim(limiti);
        set(gca, 'XTick', limiti, 'XTickLabel', [], ...
                 'YTick', limiti, 'YTickLabel', {'0.85', '1.08'}, ...
                 'TickLength', [0 0], ...
                 'FontSize', fsTick, 'LineWidth', 1, 'Layer', 'top');

        ylabel(areaY, 'FontSize', fsAsse, 'FontWeight', 'bold');
        xlabel(areaX, 'FontSize', fsAsse, 'FontWeight', 'bold');
        if p == 2, title(statoLabel); end
    end
    %% --- TABELLA STATISTICHE AL k DI RIFERIMENTO ---
    tabella_statistiche = stats_table(tabella_soggetti, statoLabel);
    fprintf('\n--- TABELLA RIASSUNTIVA TEST STATISTICI (k = %d, Wilcoxon signed-rank) ---\n', k_stat);
    disp(tabella_statistiche);

    %% --- STESSE STATISTICHE PER TUTTI I k  ---
    % 6 k x 3 confronti = 18 test, SENZA correzione per test multipli.
    % Z > 0 significa Area1 > Area2.
    tab_tutti = table();
    for ii = 1:numel(ks)
        Sk = stats_table(build_subject_table(T, ii), statoLabel);
        Sk.k = repmat(ks(ii), height(Sk), 1);
        tab_tutti = [tab_tutti; Sk]; %#ok<AGROW>
    end
    fprintf('\n--- STATISTICHE PER OGNI k (esplorativo, nessuna correzione per test multipli) ---\n');
    disp(tab_tutti(:, {'k','Confronto','n','Z','p_value','effect_size_r','mediana_diff'}));


    %% --- CURVE MEDIE NCD(k) PER AREA CON BARRE D'ERRORE ---
    figure('Color', 'w', 'Position', [100, 200, 950, 320]);
    hold on;
    plot_handles = []; labels_ok = {};
    for i = 1:numel(areas)
        Ta = T(strcmp(T.Area, areas{i}), :);
        if isempty(Ta), continue; end
        Hstack = cell2mat(Ta.NCD_vs_K);
        m  = mean(Hstack, 1, 'omitnan');
        sd = std(Hstack, 0, 1, 'omitnan');
        switch areas{i}
            case 'A1', col = color_A1; ls = ':o';
            case 'M1', col = color_M1; ls = '-s';
            case 'S1', col = color_S1; ls = '--d';
        end
        h = plot(ks, m, ls, 'Color', col, 'LineWidth', 2.5);
        errorbar(ks, m, sd, 'Color', col, 'LineWidth', 1.1, 'CapSize', 4, ...
            'LineStyle', 'none', 'HandleVisibility', 'off');
        plot_handles = [plot_handles, h]; %#ok<AGROW>
        labels_ok = [labels_ok, areas{i}]; %#ok<AGROW>
    end

    title(statoLabel, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('k',   'FontSize', 12, 'FontWeight', 'bold');
    ylabel('NCD', 'FontSize', 12, 'FontWeight', 'bold');
    grid on; box on;

    xlim([-2 70]);                 % margine a sinistra di k=2 e a destra di k=67
    xticks(ks);
    ylim([0.88 1.08]);
    yticks([0.9 0.95 1 1.05]);
    yticklabels({'0.9', '0.95', '1', '1.05'});

    set(gca, 'FontSize', 11, 'LineWidth', 1.1);
    % legend(plot_handles, labels_ok, 'Location', 'best');   % nella figura di riferimento non c'è


    %% --- TABELLA PER CANALE AL k DI RIFERIMENTO ---
    T_kstat = T;
    T_kstat.NCD_vs_K = cellfun(@(x) x(idx_k), T.NCD_vs_K);
    T_kstat.K_range  = cellfun(@(x) x(idx_k), T.K_range);

    %% --- MEDIANA E IQR PER AREA, SUI CANALI (k = k_stat) ---
    areas_stat = {'A1', 'M1', 'S1'};
    n_attesi = expectedCounts(st,:);

    riepilogo = table();
    for i = 1:numel(areas_stat)
        nm   = areas_stat{i};
        vals = T_kstat.NCD_vs_K(strcmp(T_kstat.Area, nm));
        vals = vals(~isnan(vals));

        if numel(vals) ~= n_attesi(i)
            warning('%s: trovati %d canali validi, attesi %d.', nm, numel(vals), n_attesi(i));
        end

        if isempty(vals), continue; end
        q = prctile(vals, [25 50 75]);
        riga = table({nm}, numel(vals), q(2), q(1), q(3), iqr(vals), ...
            'VariableNames', {'Area','n_canali','Mediana','Q1','Q3','IQR'});
        riepilogo = [riepilogo; riga]; %#ok<AGROW>
    end

    fprintf('\n--- MEDIANA E IQR PER AREA (k = %d, sui canali) ---\n', k_stat);
    disp(riepilogo);

    % stampa nel formato "mediana [Q1 - Q3]" per il testo/tabella del paper
    for i = 1:height(riepilogo)
        fprintf('%s (n=%d): %.3f [IQR %.3f - %.3f, ampiezza %.3f]\n', ...
            riepilogo.Area{i}, riepilogo.n_canali(i), riepilogo.Mediana(i), ...
            riepilogo.Q1(i), riepilogo.Q3(i), riepilogo.IQR(i));
    end
    %% --- UNA FIGURA PER SOGGETTO: NCD vs k, una serie per area ---
    cartellaFig = fullfile('figure_soggetti', statoLabel);
    ordineAree   = {'M1', 'S1', 'A1'};   % ordine della legenda come nella tua figura
    colArea      = struct('M1', color_M1, 'S1', color_S1, 'A1', color_A1);

    % limiti y comuni a tutti i soggetti, per poter confrontare le figure
    tuttiNCD = cell2mat(T.NCD_vs_K);
    pad      = 0.02 * (max(tuttiNCD(:), [], 'omitnan') - min(tuttiNCD(:), [], 'omitnan'));
    if pad == 0, pad = 0.01; end
    ylimComuni = [min(tuttiNCD(:), [], 'omitnan') - pad, max(tuttiNCD(:), [], 'omitnan') + pad];

    if salvaFigure && ~exist(cartellaFig, 'dir'), mkdir(cartellaFig); end

    soggetti = tabella_soggetti.SoggettoID;
    for i = 1:numel(soggetti)
        sid = soggetti{i};

        fig = figure('Color', 'w', 'Position', [100 100 560 450], ...
                     'Name', ['Subject ' sid]);
        hold on; box on; grid on;
        handles = []; etichette = {};

        for a = 1:numel(ordineAree)
            nm  = ordineAree{a};
            idx = strcmp(T.SoggettoID, sid) & strcmp(T.Area, nm);
            if ~any(idx), continue; end           % area non presente per questo soggetto

            M = cell2mat(T.NCD_vs_K(idx));   % (n canali) x (n k)
            y = mean(M, 1, 'omitnan');                   % media tra canali della stessa area

            h = plot(ks, y, 'o', 'MarkerSize', 8, 'LineStyle', 'none', ...
                     'MarkerFaceColor', colArea.(nm), 'MarkerEdgeColor', colArea.(nm));
            handles = [handles, h];          %#ok<AGROW>
            etichette = [etichette, nm];     %#ok<AGROW>
        end

        xticks(ks);
        xlim([0, max(ks) + 4]);
        ylim(ylimComuni);
        xlabel('k', 'FontWeight', 'bold');
        ylabel('NCD', 'FontWeight', 'bold');
        title(sprintf('Subject %s - %s - PPMd-H', sid, statoLabel), 'FontWeight', 'bold');
        legend(handles, etichette, 'Location', 'northeast');
        set(gca, 'FontSize', 11, 'LineWidth', 1.1);
        ax = gca; ax.XAxis.FontWeight = 'bold'; ax.YAxis.FontWeight = 'bold';

        if salvaFigure
            exportgraphics(fig, fullfile(cartellaFig, sprintf('NCD_subject_%s_%s.png', statoLabel, sid)), ...
                           'Resolution', 200);
        end
    end

    subjectTables.(statoLabel) = tabella_soggetti;
    stats_sleep = [stats_sleep; tabella_statistiche]; %#ok<AGROW>
    tab_tutti_sleep = [tab_tutti_sleep; tab_tutti]; %#ok<AGROW>
    if ~isempty(riepilogo)
        riepilogo.Stato = repmat(string(statoLabel), height(riepilogo), 1);
        riepilogo_sleep = [riepilogo_sleep; riepilogo]; %#ok<AGROW>
    end
end

fprintf('\n--- STATISTICHE TRA AREE (tutti gli stati) ---\n');
disp(stats_sleep);
fprintf('\n--- STATISTICHE ESPLORATIVE PER TUTTI I k (senza correzione) ---\n');
disp(tab_tutti_sleep);
fprintf('\n--- MEDIANA E IQR PER STATO E AREA ---\n');
disp(riepilogo_sleep);

%% --- CONFRONTI TRA FASI DEL SONNO (analoghi al confronto Awake > REM) ---
% Test appaiati a una coda, mantenuti dal precedente run del sonno.
% Le analisi tra aree, invece, sono a due code come nella veglia.
stageStats = table();
stageComparisons = {'REM','N2'; 'N2','N3'};
for a = 1:numel(areas)
    nm = areas{a};
    for c = 1:size(stageComparisons,1)
        stateHigh = stageComparisons{c,1};
        stateLow = stageComparisons{c,2};
        tabHigh = subjectTables.(stateHigh);
        tabLow = subjectTables.(stateLow);
        [~,ia,ib] = intersect(tabHigh.SoggettoID, tabLow.SoggettoID, 'stable');
        X = tabHigh.(nm)(ia);
        Y = tabLow.(nm)(ib);
        mask = ~isnan(X) & ~isnan(Y);
        X = X(mask); Y = Y(mask);
        n = numel(X);
        p = NaN; Z = NaN; md = NaN;
        if n >= 3
            [p,~,statW] = signrank(X, Y, 'tail','right');
            if isfield(statW,'zval'), Z = statW.zval; end
            md = median(X-Y);
        end
        stageStats = [stageStats; table(string(nm), ...
            string([stateHigh ' > ' stateLow]), n, p, Z, md, ...
            'VariableNames',{'Area','Comparison','n','p_value','Z','median_difference'})]; %#ok<AGROW>
    end
end
fprintf('\n--- WILCOXON TRA FASI DEL SONNO (k = %d, una coda) ---\n', k_stat);
disp(stageStats);

%% ===================== FUNZIONI LOCALI =====================
function tab = build_subject_table(T, idx_k)
% media per soggetto e per area del valore NCD all'indice idx_k
    soggetti_unici = unique(T.SoggettoID, 'stable');
    tab = table(soggetti_unici, 'VariableNames', {'SoggettoID'});
    tab.A1 = NaN(height(tab), 1);
    tab.M1 = NaN(height(tab), 1);
    tab.S1 = NaN(height(tab), 1);
    for i = 1:height(tab)
        sid = tab.SoggettoID{i};
        for area = {'A1','M1','S1'}
            idx = strcmp(T.SoggettoID, sid) & strcmp(T.Area, area{1});
            if any(idx)
                valori = cellfun(@(x) x(idx_k), T.NCD_vs_K(idx));
                tab.(area{1})(i) = mean(valori, 'omitnan');
            end
        end
    end
end

function tabS = stats_table(tab, statoLabel)
% Wilcoxon signed-rank (approssimato) tra aree, appaiato per soggetto
    Stato     = repmat(string(statoLabel), 3, 1);
    Confronto = ["S1_vs_A1"; "M1_vs_S1"; "M1_vs_A1"];
    Area1     = ["S1"; "M1"; "M1"];
    Area2     = ["A1"; "S1"; "A1"];
    n = zeros(3,1); W = NaN(3,1); p_value = NaN(3,1); Z = NaN(3,1);
    effect_size_r = NaN(3,1); mediana_diff = NaN(3,1);

    for c = 1:3
        X = tab.(char(Area2(c)));
        Y = tab.(char(Area1(c)));
        mask = ~isnan(X) & ~isnan(Y);
        Xv = X(mask); Yv = Y(mask);
        n(c) = numel(Xv);
        if n(c) >= 3
            [p_val, ~, st] = signrank(Yv, Xv, 'method', 'approximate', 'tail', 'both');
            p_value(c) = p_val;
            if isfield(st, 'signedrank'), W(c) = st.signedrank; end
            if isfield(st, 'zval'),       Z(c) = st.zval;       end
            if ~isnan(Z(c)), effect_size_r(c) = abs(Z(c)) / sqrt(n(c)); end
            mediana_diff(c) = median(Yv - Xv);
        end
    end
    tabS = table(Stato, Confronto, Area1, Area2, n, W, p_value, Z, effect_size_r, mediana_diff);
end
