% =========================================================================
%  edit_hippocampalLayers.m
%
%  Recorre cada subcarpeta de una carpeta madre, abre el archivo
%    <nombre_subcarpeta>_hippocampalLayers_channelinfo.mat
%  muestra los valores actuales y permite editarlos directamente
%  desde la consola (Enter = mantener valor actual).
%
%  Campos editables: oriens, pyramidal, radiatum, slm
%  Nota: NaN = canal no asignado (escribe NaN o deja vacio)
% =========================================================================

clear; clc;

FIELDS = {'oriens', 'pyramidal', 'radiatum', 'slm'};

%% -- 1. Seleccionar carpeta madre -----------------------------------------
parentFolder = uigetdir(pwd, 'Selecciona la carpeta madre');
if parentFolder == 0
    disp('Cancelado por el usuario.');
    return
end

%% -- 2. Listar subcarpetas ------------------------------------------------
entries    = dir(parentFolder);
subfolders = entries([entries.isdir] & ~ismember({entries.name}, {'.','..'}));
subfolders = subfolders(~startsWith({subfolders.name}, '.'));

if isempty(subfolders)
    errordlg('No se encontraron subcarpetas.', 'Error');
    return
end

fprintf('\n==========================================================\n');
fprintf('  Carpeta madre : %s\n', parentFolder);
fprintf('  Subcarpetas   : %d\n', numel(subfolders));
fprintf('==========================================================\n');
fprintf('  INTRO = mantener valor actual  |  s = saltar subcarpeta\n');
fprintf('==========================================================\n\n');

%% -- 3. Iterar subcarpetas ------------------------------------------------
for i = 1:numel(subfolders)

    subName = subfolders(i).name;
    subPath = fullfile(parentFolder, subName);
    matFile = fullfile(subPath, [subName '.hippocampalLayers.channelinfo.mat']);

    fprintf('----------------------------------------------------------\n');
    fprintf('[%d/%d]  %s\n', i, numel(subfolders), subName);

    % -- 3a. Comprobar que el archivo existe --------------------------------
    if ~isfile(matFile)
        fprintf('  AVISO: Archivo no encontrado, se omite.\n\n');
        continue
    end

    % -- 3b. Cargar --------------------------------------------------------
    S  = load(matFile, 'hippocampalLayers');
    hl = S.hippocampalLayers;

    % -- 3c. Mostrar tabla de valores actuales -----------------------------
    fprintf('\n  %-14s  %-14s  %-14s\n', 'Campo', 'layers{1,1}', 'bestShankLayers');
    fprintf('  %-14s  %-14s  %-14s\n', '--------------', '--------------', '--------------');
    for f = 1:numel(FIELDS)
        fn = FIELDS{f};
        fprintf('  %-14s  %-14s  %-14s\n', fn, ...
            formatVal(hl.layers{1,1}.(fn)), ...
            formatVal(hl.bestShankLayers.(fn)));
    end
    fprintf('\n');

    % -- 3d. Preguntar si saltar -------------------------------------------
    resp = strtrim(input('  Editar esta sesion? [INTRO=si / s=saltar]: ', 's'));
    if strcmpi(resp, 's')
        fprintf('  -> Saltada.\n\n');
        continue
    end
    fprintf('\n  Escribe el nuevo valor o pulsa INTRO para mantener el actual.\n\n');

    % -- 3e. Edicion campo a campo en consola ------------------------------
    changed = false;
    for f = 1:numel(FIELDS)
        fn = FIELDS{f};

        % layers{1,1}
        prompt_l = sprintf('  layers{1,1}.%-12s [%s]: ', fn, formatVal(hl.layers{1,1}.(fn)));
        raw_l    = strtrim(input(prompt_l, 's'));
        if ~isempty(raw_l)
            new_l = parseInput(raw_l);
            if ~isequal_nan(new_l, hl.layers{1,1}.(fn))
                fprintf('    => %s  ->  %s\n', formatVal(hl.layers{1,1}.(fn)), formatVal(new_l));
                hl.layers{1,1}.(fn) = new_l;
                changed = true;
            end
        end

        % bestShankLayers
        prompt_b = sprintf('  bestShankLayers.%-8s [%s]: ', fn, formatVal(hl.bestShankLayers.(fn)));
        raw_b    = strtrim(input(prompt_b, 's'));
        if ~isempty(raw_b)
            new_b = parseInput(raw_b);
            if ~isequal_nan(new_b, hl.bestShankLayers.(fn))
                fprintf('    => %s  ->  %s\n', formatVal(hl.bestShankLayers.(fn)), formatVal(new_b));
                hl.bestShankLayers.(fn) = new_b;
                changed = true;
            end
        end

        fprintf('\n');
    end

    % -- 3f. Guardar -------------------------------------------------------
    S.hippocampalLayers = hl;
    save(matFile, '-struct', 'S', '-append');

    if changed
        fprintf('  OK  Guardado con cambios: %s\n\n', matFile);
    else
        fprintf('  OK  Guardado sin cambios.\n\n');
    end

end

fprintf('==========================================================\n');
fprintf('  Proceso completado.\n');
fprintf('==========================================================\n');


%% =========================================================================
%  Funciones auxiliares
%% =========================================================================

function s = formatVal(v)
    if isempty(v) || (isscalar(v) && isnan(v))
        s = 'NaN';
    elseif isscalar(v)
        s = num2str(v);
    else
        s = mat2str(v(:)', 4);
    end
end

function v = parseInput(str)
    str = strtrim(str);
    if isempty(str) || strcmpi(str, 'nan')
        v = NaN;
    else
        v = str2double(str);
        if isnan(v)
            warning('edit_hippocampalLayers:badInput', ...
                'No se pudo convertir "%s" a numero; se asigna NaN.', str);
        end
    end
end

function tf = isequal_nan(a, b)
    if (isscalar(a) && isnan(a)) && (isscalar(b) && isnan(b))
        tf = true;
    else
        tf = isequal(a, b);
    end
end