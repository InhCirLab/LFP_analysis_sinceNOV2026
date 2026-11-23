function [groupeddata] = ps1to120_Alberto(thisgrouplabel, mainFolder)
% PS1TO120  Calcula espectros LFP 1-120 Hz por capa hippocampal.
%
%  - Acepta mainFolder como argumento opcional (si no, abre uigetdir).
%  - Campos delta_hires_indb / delta_hires_inuv por capa (pwelch, 1-3.5 Hz).
%  - Si falta alguna capa (or, pyr, rad, slm), la omite y procesa las demas
%    (esa sesion queda con NaN en la capa que falta).
%  - PARALELO: cada sesion se procesa en un worker (parfor). Si no hay pool
%    activo, parfor corre en serie sin problemas.
%  - Lee <sesion>.thetaEpochs6to12.states.mat (el refrescado con ripples nuevos).

% ── Carpeta principal ─────────────────────────────────────────────────────
if nargin < 2 || isempty(mainFolder)
    mainFolder = uigetdir();
end
if mainFolder == 0
    error('ps1to120: no se selecciono carpeta. Abortado.');
end

folderInfo = dir(mainFolder);
subfolders = folderInfo([folderInfo.isdir]);
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
lfppaths   = {subfolders.name};
nSess      = numel(lfppaths);

DELTA_RANGE = [1 3.5];
ALL_LAYERS  = {'or','pyr','rad','slm'};
layerFieldMap = struct('or','oriens', 'pyr','pyramidal', 'rad','radiatum', 'slm','slm');

% ── Procesar sesiones en paralelo ─────────────────────────────────────────
res = cell(nSess, 1);
parfor ii = 1:nSess
    res{ii} = processSession(fullfile(mainFolder, lfppaths{ii}), lfppaths{ii}, ...
                             ALL_LAYERS, layerFieldMap, DELTA_RANGE);
end

% ── Mensajes (workers no imprimen en el cliente) ──────────────────────────
for ii = 1:nSess
    if isempty(res{ii}), continue; end
    for m = 1:numel(res{ii}.msg), fprintf('  [%s] %s\n', lfppaths{ii}, res{ii}.msg{m}); end
end

% ── Ensamblar acumuladores sesion x frecuencia ────────────────────────────
okSess = find(~cellfun(@(r) isempty(r) || ~isfield(r,'f') || isempty(r.f), res));
if isempty(okSess)
    error('ps1to120: ninguna sesion produjo datos.');
end
freqs          = res{okSess(1)}.f;
freqs_delta_hr = res{okSess(1)}.f_delta;

varList = {'allrec_db','allrec_uv','loc_db','loc_uv','nlnr_db','nlnr_uv','nl_db','nl_uv','delta_db','delta_uv'};
acc = struct();
for li = 1:numel(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    for vi = 1:numel(varList)
        vname = varList{vi};
        % longitud del vector (primera sesion que tenga esta capa)
        len = [];
        for k = okSess(:)'
            if isfield(res{k}.layers, lname)
                len = numel(res{k}.layers.(lname).(vname));
                break
            end
        end
        if isempty(len), continue; end   % capa nunca disponible
        M = nan(nSess, len);
        for k = okSess(:)'
            if isfield(res{k}.layers, lname)
                M(k, :) = res{k}.layers.(lname).(vname);
            end
        end
        acc.(lname).(vname) = M;
    end
end

% ── ID TABLE ──────────────────────────────────────────────────────────────
ssid      = lfppaths.';
ssexpcond = repmat(thisgrouplabel, nSess, 1);
SSidtable = table(ssid, ssexpcond, 'VariableNames', ["ID", "ExpCond"]);

% ── EMPAQUETAR OUTPUT ─────────────────────────────────────────────────────
groupeddata.idtable           = SSidtable;
groupeddata.freqs             = freqs;
groupeddata.freqs_delta_hires = freqs_delta_hr;
groupeddata.label             = thisgrouplabel;

for li = 1:length(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    if ~isfield(acc, lname)
        groupeddata.(lname) = struct();
        continue
    end
    a = acc.(lname);
    groupeddata.(lname).allrec              = a.allrec_db;
    groupeddata.(lname).allrecinuv          = a.allrec_uv;
    groupeddata.(lname).locindb             = a.loc_db;
    groupeddata.(lname).locinuv             = a.loc_uv;
    groupeddata.(lname).nonlocnonrippleindb = a.nlnr_db;
    groupeddata.(lname).nonlocnonrippleinuv = a.nlnr_uv;
    groupeddata.(lname).nonlocindb          = a.nl_db;
    groupeddata.(lname).nonlocinuv          = a.nl_uv;
    groupeddata.(lname).delta_hires_indb    = a.delta_db;
    groupeddata.(lname).delta_hires_inuv    = a.delta_uv;
end
end


function r = processSession(sessPath, sessName, ALL_LAYERS, layerFieldMap, DELTA_RANGE)
% Procesa UNA sesion (corre dentro de un worker).
r = struct('f',[], 'f_delta',[], 'layers',struct());
r.msg = {};

cd(sessPath);

S = load(fullfile(sessPath, [sessName '.thetaEpochs6to12.states.mat']), 'thisthetaEpochs');
thetaEpochs = S.thisthetaEpochs;

% ── Canales disponibles ───────────────────────────────────────────────────
hippocampalLayers = getHippocampalLayers;
layerDef = hippocampalLayers.layers{1,1};

availableLayers = {};
chanMap = struct();
for li = 1:length(ALL_LAYERS)
    lname = ALL_LAYERS{li};
    fname = layerFieldMap.(lname);

    if ~isfield(layerDef, fname)
        r.msg{end+1} = sprintf('AVISO: campo "%s" inexistente, capa "%s" omitida.', fname, lname);
        continue
    end
    chNum = layerDef.(fname);
    if isempty(chNum) || ~isnumeric(chNum) || any(isnan(chNum)) || any(chNum <= 0) || any(chNum ~= floor(chNum))
        r.msg{end+1} = sprintf('AVISO: canal invalido (%s = %s) para capa "%s", omitida.', fname, mat2str(chNum), lname);
        continue
    end
    availableLayers{end+1} = lname; %#ok<AGROW>
    chanMap.(lname) = chNum;
end

if isempty(availableLayers)
    r.msg{end+1} = 'AVISO: ninguna capa disponible, sesion omitida.';
    return
end

% ── LFP de capas disponibles ──────────────────────────────────────────────
lfpData = struct();
for li = 1:length(availableLayers)
    lname = availableLayers{li};
    lfpData.(lname) = getLFP(chanMap.(lname), 'noPrompts', true);
end
Fs = lfpData.(availableLayers{1}).samplingRate;

params.Fs     = Fs;
params.fpass  = [1 120];
params.tapers = [3 5];
params.pad    = 1;

thisIDints = thetaEpochs.thetaRun.newintervals;

% ── NON-LOCOMOTION NON-RIPPLES ────────────────────────────────────────────
nlnr_int      = thetaEpochs.intervals_nonTheta_nonRipples;
dif_nlnr      = nlnr_int(:,2) - nlnr_int(:,1);
ntnr_int_keep = round(nlnr_int(dif_nlnr >= 1, :));

% ── NON-LOCOMOTION (INCLUYENDO RIPPLES) ───────────────────────────────────
nl_int      = thetaEpochs.intervals_nonTheta;
dif_nl      = nl_int(:,2) - nl_int(:,1);
nl_int_keep = round(nl_int(dif_nl >= 1, :));

% ── pwelch delta alta resolucion ──────────────────────────────────────────
nfft_hires = 2^nextpow2(Fs / 0.05);
noverlap   = floor(nfft_hires * 0.5);
win_hires  = hann(nfft_hires);
didx_hr    = [];

for li = 1:length(availableLayers)
    lname = availableLayers{li};
    lfp   = lfpData.(lname);

    % Espectrograma multitaper
    [Sx, ~, f] = mtspecgramc_fast(single(lfp.data), [2 1], params);
    Sx(Sx == 0) = NaN;
    Sdb = pow2db(Sx);

    if isempty(r.f), r.f = f; end

    % Correccion de limites para intervalos non-loc
    checkend = size(Sx, 1);
    if ntnr_int_keep(end) > checkend
        ntnr_int_keep = ntnr_int_keep(1:end-1, :);
    end
    if ntnr_int_keep(1,1) == 0
        ntnr_int_keep(1,1) = 1;
    end
    if nl_int_keep(end) > checkend
        nl_int_keep = nl_int_keep(1:end-1, :);
    end
    if nl_int_keep(1,1) == 0
        nl_int_keep(1,1) = 1;
    end

    L = struct();
    L.allrec_db = median(Sdb, 1, 'omitnan');
    L.allrec_uv = median(Sx,  1, 'omitnan');
    L.loc_db    = computeLocalizedMedian(Sdb, thisIDints);
    L.loc_uv    = computeLocalizedMedian(Sx,  thisIDints);
    L.nlnr_db   = computeLocalizedMedian(Sdb, ntnr_int_keep);
    L.nlnr_uv   = computeLocalizedMedian(Sx,  ntnr_int_keep);
    L.nl_db     = computeLocalizedMedian(Sdb, nl_int_keep);
    L.nl_uv     = computeLocalizedMedian(Sx,  nl_int_keep);

    % Delta alta resolucion (pwelch)
    sig = double(lfp.data(:,1));
    [pxx, f_hr] = pwelch(sig, win_hires, noverlap, nfft_hires, Fs);
    if isempty(didx_hr)
        didx_hr = f_hr >= DELTA_RANGE(1) & f_hr <= DELTA_RANGE(2);
        if isempty(r.f_delta), r.f_delta = f_hr(didx_hr); end
    end
    pxx_delta  = pxx(didx_hr);
    L.delta_db = pow2db(pxx_delta(:))';
    L.delta_uv = pxx_delta(:)';

    % filas fila-vector (1 x nFreq)
    fn = fieldnames(L);
    for k = 1:numel(fn), L.(fn{k}) = L.(fn{k})(:)'; end

    r.layers.(lname) = L;
end
r.msg{end+1} = 'OK';
end

% function [groupeddata] = ps1to120_Alberto(thisgrouplabel, mainFolder)
% % PS1TO120  Calcula espectros LFP 1-120 Hz por capa hippocampal.
% %
% %  CAMBIOS:
% %   1. Acepta mainFolder como argumento opcional.
% %   2. Añade campos delta_hires_indb / delta_hires_inuv por capa.
% %   3. Si falta alguna capa (or, pyr, rad, slm), la omite y procesa las demás.
% 
% % ── Carpeta principal ─────────────────────────────────────────────────────
% if nargin < 2 || isempty(mainFolder)
%     mainFolder = uigetdir();
% end
% if mainFolder == 0
%     error('ps1to120: no se seleccionó carpeta. Abortado.');
% end
% 
% folderInfo = dir(mainFolder);
% subfolders = folderInfo([folderInfo.isdir]);
% subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
% lfppaths   = {subfolders.name};
% 
% % ── Parámetros delta alta resolución ─────────────────────────────────────
% DELTA_RANGE = [1 3.5];
% freqs_delta_hr = [];   % se rellenará en la primera sesión válida
% 
% % ── Acumuladores por capa (se crean dinámicamente al primer uso) ──────────
% ALL_LAYERS = {'or','pyr','rad','slm'};
% 
% % Estructura donde guardaremos los acumuladores sesión×frecuencia
% acc = struct();   % acc.(layer).(variable)(ii,:)
% 
% for ii = 1:length(lfppaths)
% 
%     thisfoldername = strcat(mainFolder, "\", lfppaths(ii));
%     cd(thisfoldername)
%     thetafilename  = strcat(thisfoldername, "\", lfppaths(ii), ".thetaEpochs6to12.states.mat");
%     load(thetafilename);
%     thetaEpochs = thisthetaEpochs;
% 
%     % ── Obtener canales disponibles ───────────────────────────────────────
%     hippocampalLayers = getHippocampalLayers;
%     layerDef = hippocampalLayers.layers{1,1};
% 
%     % Mapa nombre → campo en layerDef
%     layerFieldMap = struct('or','oriens', 'pyr','pyramidal', 'rad','radiatum', 'slm','slm');
% 
%     availableLayers = {};
%     chanMap = struct();
%     for li = 1:length(ALL_LAYERS)
%         lname = ALL_LAYERS{li};
%         fname = layerFieldMap.(lname);
% 
%         % 1) ¿Existe el campo?
%         if ~isfield(layerDef, fname)
%             warning('ps1to120: sesión %s – campo "%s" inexistente, capa "%s" omitida.', ...
%                 lfppaths{ii}, fname, lname);
%             continue
%         end
% 
%         chNum = layerDef.(fname);
% 
%         % 2) ¿Es un entero positivo válido?
%         if isempty(chNum) || ~isnumeric(chNum) || any(chNum <= 0) || any(chNum ~= floor(chNum))
%             warning('ps1to120: sesión %s – canal inválido (%s = %s) para capa "%s", omitida.', ...
%                 lfppaths{ii}, fname, mat2str(chNum), lname);
%             continue
%         end
% 
%         availableLayers{end+1} = lname; %#ok<AGROW>
%         chanMap.(lname) = chNum;
%     end
% 
%     if isempty(availableLayers)
%         warning('ps1to120: sesión %s – ninguna capa disponible, sesión omitida.', lfppaths{ii});
%         continue
%     end
% 
%     % ── Cargar LFP solo de capas disponibles ─────────────────────────────
%     lfpData = struct();
%     for li = 1:length(availableLayers)
%         lname = availableLayers{li};
%         lfpData.(lname) = getLFP(chanMap.(lname), 'noPrompts', true);
%     end
% 
%     Fs = lfpData.(availableLayers{1}).samplingRate;
% 
%     params.Fs     = Fs;
%     params.fpass  = [1 120];
%     params.tapers = [3 5];
%     params.pad    = 1;
% 
%     thisIDints = thetaEpochs.thetaRun.newintervals;
% 
%     % ── NON-LOCOMOTION NON-RIPPLES ────────────────────────────────────────
%     nlnr_int      = thetaEpochs.intervals_nonTheta_nonRipples;
%     dif_nlnr      = nlnr_int(:,2) - nlnr_int(:,1);
%     ntnr_int_keep = round(nlnr_int(dif_nlnr >= 1, :));
% 
%     % ── NON-LOCOMOTION (INCLUDING RIPPLES) ────────────────────────────────
%     nl_int      = thisthetaEpochs.intervals_nonTheta;
%     dif_nl      = nl_int(:,2) - nl_int(:,1);
%     nl_int_keep = round(nl_int(dif_nl >= 1, :));
% 
%     % ── Parámetros pwelch delta alta resolución ───────────────────────────
%     nfft_hires = 2^nextpow2(Fs / 0.05);
%     noverlap   = floor(nfft_hires * 0.5);
%     win_hires  = hann(nfft_hires);
%     didx_hr    = [];   % se calcula en el primer canal de la primera sesión
% 
%     % ── Procesar cada capa disponible ────────────────────────────────────
%     for li = 1:length(availableLayers)
%         lname = availableLayers{li};
%         lfp   = lfpData.(lname);
% 
%         % Espectrograma multitaper
%         [S, ~, f] = mtspecgramc_fast(single(lfp.data), [2 1], params);
%         S(S == 0) = NaN;
%         Sdb       = pow2db(S);
% 
%         % Guardar vector de frecuencias la primera vez
%         if ~isfield(acc, 'freqs')
%             acc.freqs = f;
%         end
% 
%         % Corrección de límites para intervalos non-loc
%         checkend = size(S, 1);
%         if ntnr_int_keep(end) > checkend
%             ntnr_int_keep = ntnr_int_keep(1:end-1, :);
%         end
%         if ntnr_int_keep(1,1) == 0
%             ntnr_int_keep(1,1) = 1;
%         end
%         if nl_int_keep(end) > checkend
%             nl_int_keep = nl_int_keep(1:end-1, :);
%         end
%         if nl_int_keep(1,1) == 0
%             nl_int_keep(1,1) = 1;
%         end
% 
%         % Acumuladores estándar
%         acc.(lname).allrec_db(ii, :)     = median(Sdb, 1, 'omitnan');
%         acc.(lname).allrec_uv(ii, :)     = median(S,   1, 'omitnan');
%         acc.(lname).loc_db(ii, :)        = computeLocalizedMedian(Sdb, thisIDints);
%         acc.(lname).loc_uv(ii, :)        = computeLocalizedMedian(S,   thisIDints);
%         acc.(lname).nlnr_db(ii, :)       = computeLocalizedMedian(Sdb, ntnr_int_keep);
%         acc.(lname).nlnr_uv(ii, :)       = computeLocalizedMedian(S,   ntnr_int_keep);
%         acc.(lname).nl_db(ii, :)         = computeLocalizedMedian(Sdb, nl_int_keep);
%         acc.(lname).nl_uv(ii, :)         = computeLocalizedMedian(S,   nl_int_keep);
% 
%         % Delta alta resolución (pwelch)
%         sig = double(lfp.data(:,1));
%         [pxx, f_hr] = pwelch(sig, win_hires, noverlap, nfft_hires, Fs);
% 
%         if isempty(didx_hr)
%             didx_hr        = f_hr >= DELTA_RANGE(1) & f_hr <= DELTA_RANGE(2);
%             freqs_delta_hr = f_hr(didx_hr);
%         end
% 
%         pxx_delta = pxx(didx_hr);
%         acc.(lname).delta_db(ii, :) = pow2db(pxx_delta(:))';
%         acc.(lname).delta_uv(ii, :) = pxx_delta(:)';
%     end
% end
% 
% % ── ID TABLE ──────────────────────────────────────────────────────────────
% ssid      = lfppaths.';
% ssexpcond = repmat(thisgrouplabel, length(lfppaths), 1);
% SSidtable = table(ssid, ssexpcond, 'VariableNames', ["ID", "ExpCond"]);
% 
% % ── EMPAQUETAR OUTPUT ─────────────────────────────────────────────────────
% groupeddata.idtable           = SSidtable;
% groupeddata.freqs             = acc.freqs;
% groupeddata.freqs_delta_hires = freqs_delta_hr;
% groupeddata.label             = thisgrouplabel;
% 
% for li = 1:length(ALL_LAYERS)
%     lname = ALL_LAYERS{li};
%     if ~isfield(acc, lname)
%         % Capa nunca disponible en ninguna sesión: campo vacío
%         groupeddata.(lname) = struct();
%         continue
%     end
%     a = acc.(lname);
%     groupeddata.(lname).allrec              = a.allrec_db;
%     groupeddata.(lname).allrecinuv          = a.allrec_uv;
%     groupeddata.(lname).locindb             = a.loc_db;
%     groupeddata.(lname).locinuv             = a.loc_uv;
%     groupeddata.(lname).nonlocnonrippleindb = a.nlnr_db;
%     groupeddata.(lname).nonlocnonrippleinuv = a.nlnr_uv;
%     groupeddata.(lname).nonlocindb          = a.nl_db;
%     groupeddata.(lname).nonlocinuv          = a.nl_uv;
%     groupeddata.(lname).delta_hires_indb    = a.delta_db;
%     groupeddata.(lname).delta_hires_inuv    = a.delta_uv;
% end
% 
% clearvars -except groupeddata
% end
% % function [groupeddata] = ps1to120_Alberto(thisgrouplabel, mainFolder)
% % % PS1TO120  Calcula espectros LFP 1-120 Hz por capa hippocampal.
% % %
% % %  CAMBIOS respecto a la versión original:
% % %   1. Acepta mainFolder como argumento opcional. Si no se pasa (o es ''),
% % %      abre el diálogo uigetdir como antes. Esto evita el doble diálogo
% % %      cuando el script llamador ya tiene la carpeta elegida.
% % %   2. Añade campos delta_hires_indb / delta_hires_inuv por capa:
% % %      vector PSD de alta resolución en [1-3.5 Hz] calculado con pwelch
% % %      sobre el LFP completo de cada sesión (~50 puntos en delta).
% % %      freqs_delta_hires contiene el vector de frecuencias correspondiente.
% % 
% % % ── Carpeta principal ─────────────────────────────────────────────────────
% % if nargin < 2 || isempty(mainFolder)
% %     mainFolder = uigetdir();
% % end
% % if mainFolder == 0
% %     error('ps1to120: no se seleccionó carpeta. Abortado.');
% % end
% % 
% % folderInfo = dir(mainFolder);
% % subfolders = folderInfo([folderInfo.isdir]);
% % subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'}));
% % lfppaths   = {subfolders.name};
% % 
% % % ── Parámetros delta alta resolución ─────────────────────────────────────
% % DELTA_RANGE  = [1 3.5];
% % DELTA_NFFT   = 2^17;      % 131072 → df ≈ 0.0095 Hz @ Fs=1250, ~260 pts en delta
% % %  Se ajustará automáticamente a Fs real del primer registro si difiere.
% % 
% % for ii = 1:length(lfppaths)
% % 
% %     thisfoldername = strcat(mainFolder, "\", lfppaths(ii));
% %     cd(thisfoldername)
% %     thetafilename = strcat(thisfoldername, "\", lfppaths(ii), ".thetaEpochs6to12.states.mat");
% %     load(thetafilename);
% %     thetaEpochs   = thisthetaEpochs;
% %     hippocampalLayers = getHippocampalLayers;
% %     slmchannel = hippocampalLayers.layers{1,1}.slm;
% %     orchannel  = hippocampalLayers.layers{1,1}.oriens;
% %     pyrchannel = hippocampalLayers.layers{1,1}.pyramidal;
% %     radchannel = hippocampalLayers.layers{1,1}.radiatum;
% %     thisIDints = thetaEpochs.thetaRun.newintervals;
% % 
% %     slmlfpT = getLFP(slmchannel, 'noPrompts', true);
% %     orlfpT  = getLFP(orchannel,  'noPrompts', true);
% %     pyrlfpT = getLFP(pyrchannel, 'noPrompts', true);
% %     radlfpT = getLFP(radchannel, 'noPrompts', true);
% % 
% %     Fs = slmlfpT.samplingRate;
% % 
% %     slmparams.Fs     = Fs;
% %     slmparams.fpass  = [1 120];
% %     slmparams.tapers = [3 5];
% %     slmparams.pad    = 1;
% %     orparams  = slmparams;
% %     pyrparams = slmparams;
% %     radparams = slmparams;
% % 
% %     % ── OR ────────────────────────────────────────────────────────────────
% %     [orS, ort, orf] = mtspecgramc_fast(single(orlfpT.data), [2 1], orparams);
% %     orS(orS == 0)   = NaN;
% %     orSdb           = pow2db(orS);
% % 
% %     orallrecthisIDps_inuv(ii, :) = median(orS,   1, 'omitnan');
% %     orallrecthisIDps_indb(ii, :) = median(orSdb, 1, 'omitnan');
% % 
% %     % ── PYR ───────────────────────────────────────────────────────────────
% %     [pyrS, pyrt, pyrf] = mtspecgramc_fast(single(pyrlfpT.data), [2 1], pyrparams);
% %     pyrS(pyrS == 0)    = NaN;
% %     pyrSdb             = pow2db(pyrS);
% % 
% %     pyrallrecthisIDps_inuv(ii, :) = median(pyrS,   1, 'omitnan');
% %     pyrallrecthisIDps_indb(ii, :) = median(pyrSdb, 1, 'omitnan');
% % 
% %     % ── RAD ───────────────────────────────────────────────────────────────
% %     [radS, radt, radf] = mtspecgramc_fast(single(radlfpT.data), [2 1], radparams);
% %     radS(radS == 0)    = NaN;
% %     radSdb             = pow2db(radS);
% % 
% %     radallrecthisIDps_inuv(ii, :) = median(radS,   1, 'omitnan');
% %     radallrecthisIDps_indb(ii, :) = median(radSdb, 1, 'omitnan');
% % 
% %     % ── SLM ───────────────────────────────────────────────────────────────
% %     [slmS, slmt, slmf] = mtspecgramc_fast(single(slmlfpT.data), [2 1], slmparams);
% %     slmS(slmS == 0)    = NaN;
% %     slmSdb             = pow2db(slmS);
% % 
% %     slmallrecthisIDps_inuv(ii, :) = median(slmS,   1, 'omitnan');
% %     slmallrecthisIDps_indb(ii, :) = median(slmSdb, 1, 'omitnan');
% % 
% %     % ── DELTA ALTA RESOLUCIÓN (pwelch sobre LFP completo) ─────────────────
% %     %  nfft ajustado para dar ~50 puntos en [1-3.5 Hz]
% %     %  df_target = 2.5 Hz / 50 puntos = 0.05 Hz
% %     %  nfft_needed = Fs / df_target
% %     nfft_hires = 2^nextpow2(Fs / 0.05);   % ≥ Fs/0.05, potencia de 2
% %     noverlap   = floor(nfft_hires * 0.5);  % 50% solapamiento
% %     win_hires  = hann(nfft_hires);
% % 
% %     channels   = {orlfpT, pyrlfpT, radlfpT, slmlfpT};
% %     chNames    = {'or','pyr','rad','slm'};
% % 
% %     for ch = 1:4
% %         sig = double(channels{ch}.data(:,1));
% %         [pxx, f_hr] = pwelch(sig, win_hires, noverlap, nfft_hires, Fs);
% %         % Guardar solo el rango delta
% %         if ii == 1 && ch == 1
% %             didx_hr       = f_hr >= DELTA_RANGE(1) & f_hr <= DELTA_RANGE(2);
% %             freqs_delta_hr = f_hr(didx_hr);
% %         end
% %         pxx_delta = pxx(didx_hr);
% %         eval(sprintf('%s_delta_hires_uv(ii,:) = pxx_delta(:)'';', chNames{ch}));
% %         eval(sprintf('%s_delta_hires_db(ii,:) = pow2db(pxx_delta(:))'';', chNames{ch}));
% %     end
% % 
% %     % ── LOCOMOTION ────────────────────────────────────────────────────────
% %     orlocS_median_db(ii, :)  = computeLocalizedMedian(orSdb,  thisIDints);
% %     pyrlocS_median_db(ii, :) = computeLocalizedMedian(pyrSdb, thisIDints);
% %     radlocS_median_db(ii, :) = computeLocalizedMedian(radSdb, thisIDints);
% %     slmlocS_median_db(ii, :) = computeLocalizedMedian(slmSdb, thisIDints);
% % 
% %     orlocS_median_uv(ii, :)  = computeLocalizedMedian(orS,  thisIDints);
% %     pyrlocS_median_uv(ii, :) = computeLocalizedMedian(pyrS, thisIDints);
% %     radlocS_median_uv(ii, :) = computeLocalizedMedian(radS, thisIDints);
% %     slmlocS_median_uv(ii, :) = computeLocalizedMedian(slmS, thisIDints);
% % 
% %     % ── NON-LOCOMOTION NON-RIPPLES ────────────────────────────────────────
% %     nlnr_int       = thetaEpochs.intervals_nonTheta_nonRipples;
% %     dif_nlnr_int   = nlnr_int(:,2) - nlnr_int(:,1);
% %     ntnr_int_keep  = nlnr_int(dif_nlnr_int >= 1, :);
% %     ntnr_int_keep  = round(ntnr_int_keep);
% %     checkend       = size(slmS, 1);
% %     if ntnr_int_keep(end) > checkend
% %         ntnr_int_keep = ntnr_int_keep(1:end-1, :);
% %     end
% %     if ntnr_int_keep(1,1) == 0
% %         ntnr_int_keep(1,1) = 1;
% %     end
% % 
% %     ornlnr_median_db(ii, :)  = computeLocalizedMedian(orSdb,  ntnr_int_keep);
% %     pyrnlnr_median_db(ii, :) = computeLocalizedMedian(pyrSdb, ntnr_int_keep);
% %     radnlnr_median_db(ii, :) = computeLocalizedMedian(radSdb, ntnr_int_keep);
% %     slmnlnr_median_db(ii, :) = computeLocalizedMedian(slmSdb, ntnr_int_keep);
% % 
% %     ornlnr_median_uv(ii, :)  = computeLocalizedMedian(orS,  ntnr_int_keep);
% %     pyrnlnr_median_uv(ii, :) = computeLocalizedMedian(pyrS, ntnr_int_keep);
% %     radnlnr_median_uv(ii, :) = computeLocalizedMedian(radS, ntnr_int_keep);
% %     slmnlnr_median_uv(ii, :) = computeLocalizedMedian(slmS, ntnr_int_keep);
% % 
% %     % ── NON-LOCOMOTION (INCLUDING RIPPLES) ────────────────────────────────
% %     nl_int       = thisthetaEpochs.intervals_nonTheta;
% %     dif_nl_int   = nl_int(:,2) - nl_int(:,1);
% %     nl_int_keep  = nl_int(dif_nl_int >= 1, :);
% %     nl_int_keep  = round(nl_int_keep);
% %     if nl_int_keep(end) > checkend
% %         nl_int_keep = nl_int_keep(1:end-1, :);
% %     end
% %     if nl_int_keep(1,1) == 0
% %         nl_int_keep(1,1) = 1;
% %     end
% % 
% %     ornl_median_db(ii, :)  = computeLocalizedMedian(orSdb,  nl_int_keep);
% %     pyrnl_median_db(ii, :) = computeLocalizedMedian(pyrSdb, nl_int_keep);
% %     radnl_median_db(ii, :) = computeLocalizedMedian(radSdb, nl_int_keep);
% %     slmnl_median_db(ii, :) = computeLocalizedMedian(slmSdb, nl_int_keep);
% % 
% %     ornl_median_uv(ii, :)  = computeLocalizedMedian(orS,  nl_int_keep);
% %     pyrnl_median_uv(ii, :) = computeLocalizedMedian(pyrS, nl_int_keep);
% %     radnl_median_uv(ii, :) = computeLocalizedMedian(radS, nl_int_keep);
% %     slmnl_median_uv(ii, :) = computeLocalizedMedian(slmS, nl_int_keep);
% % 
% % end
% % 
% % freqs = orf;
% % 
% % % ── ID TABLE ──────────────────────────────────────────────────────────────
% % ssid       = lfppaths.';
% % ssexpcond  = repmat(thisgrouplabel, length(lfppaths), 1);
% % SSidtable  = table(ssid, ssexpcond, 'VariableNames', ["ID", "ExpCond"]);
% % 
% % % ── EMPAQUETAR OUTPUT ─────────────────────────────────────────────────────
% % 
% % groupeddata.idtable         = SSidtable;
% % groupeddata.freqs           = freqs;
% % groupeddata.freqs_delta_hires = freqs_delta_hr;   % vector alta res delta
% % groupeddata.label           = thisgrouplabel;
% % 
% % % OR
% % groupeddata.or.allrec              = orallrecthisIDps_indb;
% % groupeddata.or.allrecinuv          = orallrecthisIDps_inuv;
% % groupeddata.or.locindb             = orlocS_median_db;
% % groupeddata.or.locinuv             = orlocS_median_uv;
% % groupeddata.or.nonlocnonrippleindb = ornlnr_median_db;
% % groupeddata.or.nonlocnonrippleinuv = ornlnr_median_uv;
% % groupeddata.or.nonlocindb          = ornl_median_db;
% % groupeddata.or.nonlocinuv          = ornl_median_uv;
% % groupeddata.or.delta_hires_indb    = or_delta_hires_db;   % ← NUEVO alta res
% % groupeddata.or.delta_hires_inuv    = or_delta_hires_uv;   % ← NUEVO alta res
% % 
% % % PYR
% % groupeddata.pyr.allrec              = pyrallrecthisIDps_indb;
% % groupeddata.pyr.allrecinuv          = pyrallrecthisIDps_inuv;
% % groupeddata.pyr.locindb             = pyrlocS_median_db;
% % groupeddata.pyr.locinuv             = pyrlocS_median_uv;
% % groupeddata.pyr.nonlocnonrippleindb = pyrnlnr_median_db;
% % groupeddata.pyr.nonlocnonrippleinuv = pyrnlnr_median_uv;
% % groupeddata.pyr.nonlocindb          = pyrnl_median_db;
% % groupeddata.pyr.nonlocinuv          = pyrnl_median_uv;
% % groupeddata.pyr.delta_hires_indb    = pyr_delta_hires_db;
% % groupeddata.pyr.delta_hires_inuv    = pyr_delta_hires_uv;
% % 
% % % RAD
% % groupeddata.rad.allrec              = radallrecthisIDps_indb;
% % groupeddata.rad.allrecinuv          = radallrecthisIDps_inuv;
% % groupeddata.rad.locindb             = radlocS_median_db;
% % groupeddata.rad.locinuv             = radlocS_median_uv;
% % groupeddata.rad.nonlocnonrippleindb = radnlnr_median_db;
% % groupeddata.rad.nonlocnonrippleinuv = radnlnr_median_uv;
% % groupeddata.rad.nonlocindb          = radnl_median_db;
% % groupeddata.rad.nonlocinuv          = radnl_median_uv;
% % groupeddata.rad.delta_hires_indb    = rad_delta_hires_db;
% % groupeddata.rad.delta_hires_inuv    = rad_delta_hires_uv;
% % 
% % % SLM
% % groupeddata.slm.allrec              = slmallrecthisIDps_indb;
% % groupeddata.slm.allrecinuv          = slmallrecthisIDps_inuv;
% % groupeddata.slm.locindb             = slmlocS_median_db;
% % groupeddata.slm.locinuv             = slmlocS_median_uv;
% % groupeddata.slm.nonlocnonrippleindb = slmnlnr_median_db;
% % groupeddata.slm.nonlocnonrippleinuv = slmnlnr_median_uv;
% % groupeddata.slm.nonlocindb          = slmnl_median_db;
% % groupeddata.slm.nonlocinuv          = slmnl_median_uv;
% % groupeddata.slm.delta_hires_indb    = slm_delta_hires_db;
% % groupeddata.slm.delta_hires_inuv    = slm_delta_hires_uv;
% % 
% % clearvars -except groupeddata
% % end