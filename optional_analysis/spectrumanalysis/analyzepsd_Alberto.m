%% ========================================================================
%  ANÁLISIS DE PSD — Parte IV (script independiente)
%  ========================================================================
%
%  Calcula el Power Spectral Density (PSD) en todas las bandas [0.1-120 Hz]
%  por grupo experimental y capa hippocampal, en dB Y en uV^2/Hz.
%
%  FLUJO DE DIÁLOGOS (todos al inicio, antes de cargar datos):
%    1. Número de grupos
%    2. Por cada grupo: etiqueta + carpeta de datos
%    3. Carpeta OUTPUT
%    → Solo entonces llama a ps1to120() y procesa.
%
%  SALIDA (carpeta output elegida):
%    delta_psd_<grupo>_<cond>_dB.xlsx      — vector PSD en dB  [nAnim x nFreqs, hoja por capa]
%    delta_psd_<grupo>_<cond>_uV.xlsx      — vector PSD en uV^2/Hz [nAnim x nFreqs, hoja por capa]
%    delta_animal_summary_<cond>.xlsx      — peak fq + power (dB y uV) por animal y capa
%    delta_group_summary_<cond>.xlsx       — media ± SD (dB y uV) por grupo y capa
%    figures/delta/
%      delta_PSD_<cond>_dB.png/.fig        — curvas PSD grupales en dB (media ± SEM)
%      delta_PSD_<cond>_uV.png/.fig        — curvas PSD grupales en uV^2/Hz (media ± SEM)
%      delta_bars_<cond>_dB.png/.fig       — barras peak fq + mean power en dB
%      delta_bars_<cond>_uV.png/.fig       — barras peak fq + mean power en uV^2/Hz
%
%  INSTRUCCIONES:
%    1. Edita D0 si necesitas cambiar parámetros o activar/desactivar condiciones
%    2. F5 — todos los diálogos aparecen al inicio y luego corre automáticamente
%    3. ps1to120() debe estar en el path de MATLAB
% =========================================================================

clearvars; close all; clc;

%% ── D0  PARÁMETROS ───────────────────────────────────────────────────────

DELTA_RANGE = [0.1 120];                       % Hz — banda delta
LAYER_NAMES = {'or', 'pyr', 'rad', 'slm'};

% Condiciones — comenta las filas que no necesites
% Formato: {campo_dB,  campo_uV,  etiqueta_archivos}
%   Si campo_uV no existe en la estructura de ps1to120() se convierte
%   automáticamente con db2pow() desde el campo dB.
CONDITIONS = {
    'allrec',              'allrecinuv',           'AllRec';
    'locindb',             'locinuv',              'Loc';
    'nonlocnonrippleindb', 'nonlocnonrippleinuv',  'NonLocNonRip';
    'nonlocindb',          'nonlocinuv',           'NonLoc';
};

%% ── D1  TODOS LOS DIÁLOGOS AL INICIO ────────────────────────────────────
%  Etiquetas + carpetas de grupos + carpeta output, todo ANTES de ps1to120()

fprintf('=== PARTE IV: Análisis de Delta [%.1f–%.1f Hz] ===\n', DELTA_RANGE(1), DELTA_RANGE(2));
fprintf('    Configurando grupos antes de cargar datos...\n\n');

% Número de grupos
nGrpsStr = inputdlg('¿Cuántos grupos experimentales?', 'Delta — Grupos', 1, {'2'});
if isempty(nGrpsStr), error('Abortado por el usuario.'); end
nGroups = str2double(nGrpsStr{1});
assert(~isnan(nGroups) && nGroups >= 1, 'Número de grupos no válido.');

% Etiqueta + carpeta por grupo
groupLabels  = cell(nGroups, 1);
groupFolders = cell(nGroups, 1);

for g = 1:nGroups
    lbl = inputdlg(sprintf('Etiqueta del grupo %d  (ej: WT, KO, PV...):', g), ...
                   sprintf('Grupo %d — Etiqueta', g), 1, {sprintf('G%d', g)});
    if isempty(lbl), error('Abortado por el usuario.'); end
    groupLabels{g} = strtrim(lbl{1});

    fprintf('>>> Selecciona la carpeta de datos del grupo "%s":\n', groupLabels{g});
    fld = uigetdir(pwd, sprintf('Datos — %s', groupLabels{g}));
    if fld == 0, error('No se seleccionó carpeta. Abortado.'); end
    groupFolders{g} = fld;
end

% Carpeta output
fprintf('>>> Selecciona la carpeta OUTPUT:\n');
outFolder = uigetdir(groupFolders{1}, 'OUTPUT — Análisis Delta');
if outFolder == 0, error('Abortado por el usuario.'); end

figFolder = fullfile(outFolder, 'figures', 'delta');
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

% Confirmación antes de empezar a cargar
fprintf('\n── Configuración ────────────────────────────────────────────────\n');
for g = 1:nGroups
    fprintf('  Grupo %d: "%s"  →  %s\n', g, groupLabels{g}, groupFolders{g});
end
fprintf('  Output  :  %s\n', outFolder);
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('Cargando datos con ps1to120()...\n\n');

%% ── D2  CARGAR DATOS CON ps1to120() ─────────────────────────────────────
%  Solo aquí, con toda la configuración ya confirmada

groupData = cell(nGroups, 1);
origDir   = pwd;

for g = 1:nGroups
    fprintf('  >> Grupo "%s"...\n', groupLabels{g});
    cd(groupFolders{g});
    groupData{g} = ps1to120_Alberto(groupLabels{g});
    cd(origDir);
    fprintf('     %d animal(es) cargado(s).\n\n', height(groupData{g}.idtable));
end

%% ── D3  VECTOR DE FRECUENCIAS DELTA ─────────────────────────────────────

freqs    = groupData{1}.freqs;
deltaIdx = freqs >= DELTA_RANGE(1) & freqs <= DELTA_RANGE(2);
freqVec  = freqs(deltaIdx);
nFreqs   = sum(deltaIdx);
fprintf('Vector delta: %d bins  (%.4f – %.4f Hz)\n\n', nFreqs, freqVec(1), freqVec(end));

%% ── D4  BUCLE PRINCIPAL: CONDICIONES ────────────────────────────────────

for cc = 1:size(CONDITIONS, 1)

    condFieldDB = CONDITIONS{cc, 1};
    condFieldUV = CONDITIONS{cc, 2};
    condLabel   = CONDITIONS{cc, 3};

    fprintf('══ Condición: %s ══\n', condLabel);

    animalSummary = table();
    groupSummary  = table();
    allPSD_dB     = cell(nGroups, 1);
    allPSD_uV     = cell(nGroups, 1);

    for g = 1:nGroups

        gData  = groupData{g};
        gLabel = groupLabels{g};
        nAnim  = height(gData.idtable);

        subDirs = dir(groupFolders{g});
        subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name}, {'.','..'}));
        folderNames = {subDirs.name}';
        
        if numel(folderNames) == nAnim
            animalIDs = folderNames;
        else
            warning('Número de subcarpetas (%d) ≠ animales en idtable (%d). Usando índices.', ...
                    numel(folderNames), nAnim);
            animalIDs = arrayfun(@(x) sprintf('Animal_%02d',x), (1:nAnim)', 'UniformOutput', false);
        end

        % Verificar que existe el campo dB
        if ~isfield(gData.(LAYER_NAMES{1}), condFieldDB)
            fprintf('  [AVISO] Campo "%s" no encontrado en grupo "%s". Saltando.\n', condFieldDB, gLabel);
            continue;
        end

        % ── Extraer PSD delta [nAnim x nFreqs x nLayers] en dB y uV ─────
        psdDB = nan(nAnim, nFreqs, numel(LAYER_NAMES));
        psdUV = nan(nAnim, nFreqs, numel(LAYER_NAMES));

        for L = 1:numel(LAYER_NAMES)
            % dB
            matDB = gData.(LAYER_NAMES{L}).(condFieldDB);
            psdDB(:, :, L) = matDB(:, deltaIdx);

            % uV^2/Hz — campo directo o conversión desde dB
            if isfield(gData.(LAYER_NAMES{L}), condFieldUV)
                matUV = gData.(LAYER_NAMES{L}).(condFieldUV);
                psdUV(:, :, L) = matUV(:, deltaIdx);
            else
                psdUV(:, :, L) = db2pow(psdDB(:, :, L));
                if L == 1
                    fprintf('  [INFO] "%s" no existe — uV calculado via db2pow desde dB.\n', condFieldUV);
                end
            end
        end

        allPSD_dB{g} = psdDB;
        allPSD_uV{g} = psdUV;

        % ── Guardar PSD vectorial en Excel ────────────────────────────────
        colNames = arrayfun(@(f) sprintf('%.4fHz',f), freqVec, 'UniformOutput', false);
        colNames = matlab.lang.makeValidName(colNames);
        idTab    = table(animalIDs, 'VariableNames', "Animal");

        xlsxDB = fullfile(outFolder, sprintf('delta_psd_%s_%s_dB.xlsx', gLabel, condLabel));
        xlsxUV = fullfile(outFolder, sprintf('delta_psd_%s_%s_uV.xlsx', gLabel, condLabel));

        for L = 1:numel(LAYER_NAMES)
            sh = upper(LAYER_NAMES{L});
            writetable([idTab  array2table(psdDB(:,:,L), 'VariableNames', colNames)], xlsxDB, 'Sheet', sh);
            writetable([idTab  array2table(psdUV(:,:,L), 'VariableNames', colNames)], xlsxUV, 'Sheet', sh);
        end
        fprintf('  PSD guardado: delta_psd_%s_%s_dB.xlsx  /  ..._uV.xlsx\n', gLabel, condLabel);

        % ── Métricas por animal ───────────────────────────────────────────
        peakFqBlock = nan(nAnim, numel(LAYER_NAMES));
        powPeakDB   = nan(nAnim, numel(LAYER_NAMES));
        powPeakUV   = nan(nAnim, numel(LAYER_NAMES));
        meanPowDB   = nan(nAnim, numel(LAYER_NAMES));
        meanPowUV   = nan(nAnim, numel(LAYER_NAMES));

        for L = 1:numel(LAYER_NAMES)
            mDB = psdDB(:, :, L);   % [nAnim x nFreqs]
            mUV = psdUV(:, :, L);
            [~, maxIdx]        = max(mDB, [], 2);
            linIdx             = sub2ind(size(mDB), (1:nAnim)', maxIdx);
            peakFqBlock(:, L)  = freqVec(maxIdx);
            powPeakDB(:, L)    = mDB(linIdx);
            powPeakUV(:, L)    = mUV(linIdx);
            meanPowDB(:, L)    = mean(mDB, 2, 'omitnan');
            meanPowUV(:, L)    = mean(mUV, 2, 'omitnan');
        end

        % Tabla por animal
        gAnimalTab = table( ...
            repmat({gLabel},    nAnim, 1), ...
            repmat({condLabel}, nAnim, 1), ...
            animalIDs, ...
            peakFqBlock(:,1),  peakFqBlock(:,2),  peakFqBlock(:,3),  peakFqBlock(:,4), ...
            powPeakDB(:,1),    powPeakDB(:,2),    powPeakDB(:,3),    powPeakDB(:,4), ...
            powPeakUV(:,1),    powPeakUV(:,2),    powPeakUV(:,3),    powPeakUV(:,4), ...
            meanPowDB(:,1),    meanPowDB(:,2),    meanPowDB(:,3),    meanPowDB(:,4), ...
            meanPowUV(:,1),    meanPowUV(:,2),    meanPowUV(:,3),    meanPowUV(:,4), ...
            'VariableNames', ["Group","Condition","Animal", ...
                "peakFq_or","peakFq_pyr","peakFq_rad","peakFq_slm", ...
                "powPeak_dB_or","powPeak_dB_pyr","powPeak_dB_rad","powPeak_dB_slm", ...
                "powPeak_uV_or","powPeak_uV_pyr","powPeak_uV_rad","powPeak_uV_slm", ...
                "meanPow_dB_or","meanPow_dB_pyr","meanPow_dB_rad","meanPow_dB_slm", ...
                "meanPow_uV_or","meanPow_uV_pyr","meanPow_uV_rad","meanPow_uV_slm"]);
        animalSummary = [animalSummary; gAnimalTab]; %#ok<AGROW>

        % Fila de resumen grupal (media ± SD)
        rowG = buildGroupRow(gLabel, condLabel, nAnim, peakFqBlock, ...
                             meanPowDB, meanPowUV, powPeakDB, powPeakUV);
        groupSummary = [groupSummary; rowG]; %#ok<AGROW>

    end  % grupos

    % ── Guardar tablas resumen ────────────────────────────────────────────
    writetable(animalSummary, fullfile(outFolder, sprintf('delta_animal_summary_%s.xlsx', condLabel)));
    writetable(groupSummary,  fullfile(outFolder, sprintf('delta_group_summary_%s.xlsx',  condLabel)));
    fprintf('  Tablas resumen guardadas.\n');

    % ── Figuras PSD (dB y uV) ─────────────────────────────────────────────
    saveDeltaPSDFig(allPSD_dB, freqVec, groupLabels, LAYER_NAMES, condLabel, 'dB',    figFolder);
    saveDeltaPSDFig(allPSD_uV, freqVec, groupLabels, LAYER_NAMES, condLabel, 'uV2Hz', figFolder);

    % ── Figuras barras (dB y uV) ──────────────────────────────────────────
    saveDeltaBarFig(animalSummary, groupLabels, LAYER_NAMES, condLabel, 'dB', figFolder);
    saveDeltaBarFig(animalSummary, groupLabels, LAYER_NAMES, condLabel, 'uV', figFolder);

    fprintf('  Condición "%s" completada.\n\n', condLabel);

end  % condiciones

fprintf('[Parte IV] Análisis de delta completado.\nResultados en: %s\n', outFolder);

%% ========================================================================
%  FUNCIONES AUXILIARES
%% ========================================================================

function row = buildGroupRow(gLabel, condLabel, nAnim, peakFqBlock, ...
                              meanPowDB, meanPowUV, powPeakDB, powPeakUV)
%BUILDGROUPROW  Fila con media ± SD de todas las métricas para un grupo.

    layNames = {'or','pyr','rad','slm'};
    nL       = numel(layNames);
    metrics  = {peakFqBlock, meanPowDB, meanPowUV, powPeakDB, powPeakUV};
    prefixes = {'peakFq', 'meanPow_dB', 'meanPow_uV', 'powPeak_dB', 'powPeak_uV'};

    varNames = ["Group","Condition","n"];
    vals     = {{gLabel}, {condLabel}, nAnim};

    for m = 1:numel(metrics)
        for l = 1:nL
            varNames(end+1) = sprintf('mean_%s_%s', prefixes{m}, layNames{l}); %#ok<AGROW>
            varNames(end+1) = sprintf('sd_%s_%s',   prefixes{m}, layNames{l}); %#ok<AGROW>
            vals{end+1}     = mean(metrics{m}(:,l), 'omitnan');               %#ok<AGROW>
            vals{end+1}     = std( metrics{m}(:,l), 'omitnan');               %#ok<AGROW>
        end
    end

    row = table(vals{:}, 'VariableNames', varNames);
end


function saveDeltaPSDFig(allPSD_byGroup, freqVec, groupLabels, layerNames, ...
                          condLabel, unitTag, figFolder)
%SAVADELTAPSDFIG  Curvas PSD delta (media ± SEM shading) superpuestas por grupo.
%   unitTag: 'dB' o 'uV2Hz'

    if strcmp(unitTag, 'dB')
        yUnit = 'PSD (dB)';
    else
        yUnit = 'PSD (uV^2/Hz)';
    end

    nG  = numel(groupLabels);
    nL  = numel(layerNames);
    pal = lines(nG);

    fig = figure('Visible','off','Color','w','Position',[50 50 700 220*nL]);

    for L = 1:nL
        ax = subplot(nL, 1, L);
        hold(ax, 'on');

        for g = 1:nG
            psd3d = allPSD_byGroup{g};
            if isempty(psd3d), continue; end
            mat = psd3d(:, :, L);
            mu  = mean(mat, 1, 'omitnan');
            sem = std(mat, 0, 1, 'omitnan') ./ sqrt(max(sum(~isnan(mat), 1), 1));

            xF = [freqVec(:)', fliplr(freqVec(:)')];
            yF = [mu+sem,      fliplr(mu-sem)];
            fill(ax, xF, yF, pal(g,:), 'FaceAlpha', 0.20, 'EdgeColor', 'none');
            plot(ax, freqVec, mu, '-', 'Color', pal(g,:), 'LineWidth', 2.2, ...
                 'DisplayName', groupLabels{g});
        end

        set(ax, 'FontSize',9, 'TickDir','out', 'XLim',[freqVec(1) freqVec(end)]);
        xlabel(ax, 'Frequency (Hz)', 'FontSize', 9);
        ylabel(ax, sprintf('%s\n%s', upper(layerNames{L}), yUnit), 'FontSize', 9);
        grid(ax, 'on'); box(ax, 'off');

        if L == 1
            legend(ax, 'Location','best', 'FontSize',8);
            title(ax, sprintf('Delta PSD [%.1f–%.1f Hz] — %s — %s', ...
                  freqVec(1), freqVec(end), condLabel, unitTag), ...
                  'FontSize',11, 'FontWeight','bold');
        end
    end

    fname = fullfile(figFolder, sprintf('delta_PSD_%s_%s.png', condLabel, unitTag));
    exportgraphics(fig, fname, 'Resolution', 200);
    savefig(fig, strrep(fname, '.png', '.fig'));
    close(fig);
    fprintf('  Figura PSD (%s) guardada.\n', unitTag);
end


function saveDeltaBarFig(T, groupLabels, layerNames, condLabel, unitTag, figFolder)
%SAVEDELTABARFIG  Barras media+SEM+scatter: peak fq y mean power por capa.
%   unitTag: 'dB' o 'uV'

    if strcmp(unitTag, 'dB')
        powColPfx = 'meanPow_dB';
        powLbl    = 'Mean Power (dB)';
    else
        powColPfx = 'meanPow_uV';
        powLbl    = 'Mean Power (uV^2/Hz)';
    end

    nG  = numel(groupLabels);
    nL  = numel(layerNames);
    pal = lines(nG);

    colDefs = {
        'peakFq',   'Peak Frequency (Hz)';
        powColPfx,  powLbl;
    };

    fig = figure('Visible','off','Color','w','Position',[50 50 560 210*nL]);

    for L = 1:nL
        for m = 1:2
            ax      = subplot(nL, 2, (L-1)*2 + m);
            colName = sprintf('%s_%s', colDefs{m,1}, layerNames{L});
            hold(ax,'on');

            for g = 1:nG
                mask = strcmp(T.Group, groupLabels{g});
                vals = T.(colName)(mask);
                nAn  = sum(mask);
                mu   = mean(vals, 'omitnan');
                sem  = std(vals,  'omitnan') / sqrt(nAn);

                bar(ax, g, mu, 0.55, 'FaceColor',pal(g,:), 'FaceAlpha',0.60, 'EdgeColor','none');
                errorbar(ax, g, mu, sem, 'k-', 'LineWidth',1.5, 'CapSize',6);
                scatter(ax, g+(rand(nAn,1)-0.5)*0.2, vals, 28, pal(g,:), 'filled', ...
                        'MarkerFaceAlpha',0.85, 'MarkerEdgeColor','k', 'LineWidth',0.4);
            end

            set(ax, 'XTick',1:nG, 'XTickLabel',groupLabels, 'FontSize',9, ...
                    'TickDir','out', 'XLim',[0.4, nG+0.6]);
            grid(ax,'on'); box(ax,'off');
            if L == 1, title(ax, colDefs{m,2}, 'FontSize',10,'FontWeight','bold'); end
            if m == 1, ylabel(ax, upper(layerNames{L}), 'FontSize',9); end
        end
    end

    sgtitle(sprintf('Delta Analysis — %s — %s', condLabel, unitTag), ...
            'FontSize',12,'FontWeight','bold');
    fname = fullfile(figFolder, sprintf('delta_bars_%s_%s.png', condLabel, unitTag));
    exportgraphics(fig, fname, 'Resolution', 200);
    savefig(fig, strrep(fname, '.png', '.fig'));
    close(fig);
    fprintf('  Figura barras (%s) guardada.\n', unitTag);
end