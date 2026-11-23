function [Tout, report] = flagOutliersPerMeasure(Tin, varargin)
% flagOutliersPerMeasure  Detect outliers per measure and replace with NaN.
%
% Rows are NEVER removed. Only the specific (row, measure) combination is
% set to NaN if it exceeds the robust z-score threshold.
%
% USAGE:
%   [Tout, report] = flagOutliersPerMeasure(Tin)
%   [Tout, report] = flagOutliersPerMeasure(Tin, 'thresh', 4.0, ...)
%
% INPUTS:
%   Tin         - Table. Each row = one unit. Must contain numeric measure
%                 columns plus grouping columns (subject, expcond, ID, etc.)
%
% OPTIONAL NAME-VALUE PARAMETERS:
%   'measureCols'  - Numeric indices or cell array of column names to check.
%                    Default: all numeric columns.
%   'logVars'      - Cell array of measure names to log10-transform before
%                    computing MAD. Only include strictly positive measures
%                    (e.g. firingRate). Default: {} (no log transform).
%   'thresh'       - Robust z-score threshold. Default: 4.0.
%   'minN_rec'     - Min samples in a recording to use recording-level MAD.
%                    Default: 20.
%   'minN_sub'     - Min samples at subject level. Default: 20.
%   'minN_group'   - Min samples at group level. Default: 20.
%   'recCol'       - Column name identifying recording/session. Default: 'ID'.
%   'subCol'       - Column name identifying subject. Default: 'subject'.
%   'groupCol'     - Column name identifying group. Default: 'expcond'.
%   'verbose'      - Print per-measure summary. Default: true.
%
% OUTPUTS:
%   Tout    - Same table as Tin with outlier values replaced by NaN.
%             No rows are removed. All non-measure columns are unchanged.
%   report  - Struct with fields:
%               .perMeasure : struct array, one per measure, with fields:
%                               .name       : measure name
%                               .nFlagged   : number of values set to NaN
%                               .flaggedRows: row indices flagged
%                               .level      : which fallback level was used
%                                             per flagged value
%               .totalFlagged : total (row,measure) pairs set to NaN
%               .params       : parameters used
%
% FALLBACK CASCADE (per recording, per variable):
%   If a recording has >= minN_rec valid samples  → use recording-level MAD
%   Else if subject has >= minN_sub valid samples → use subject-level MAD
%   Else if group has >= minN_group valid samples → use group-level MAD
%   Else                                          → use global MAD
%
% NOTES:
%   - Negative and zero values are excluded from log-transformed MAD
%     computation but are NOT automatically flagged. Their robust z is
%     computed on the original (untransformed) scale instead.
%   - Variables with a MAD of 0 fall back to std; if std is also 0 (all
%     identical values) no flagging is performed for that variable/level.
%   - The function does not modify non-numeric columns.

%% ---- Parse inputs ------------------------------------------------------

p = inputParser;
addRequired(p,  'Tin');
addParameter(p, 'measureCols',  [],     @(x) isnumeric(x)||iscell(x));
addParameter(p, 'logVars',      {},     @iscell);
addParameter(p, 'thresh',       4.0,   @(x) isnumeric(x) && x > 0);
addParameter(p, 'minN_rec',     20,    @(x) isnumeric(x) && x > 0);
addParameter(p, 'minN_sub',     20,    @(x) isnumeric(x) && x > 0);
addParameter(p, 'minN_group',   20,    @(x) isnumeric(x) && x > 0);
addParameter(p, 'recCol',       'ID',  @ischar);
addParameter(p, 'subCol',       'subject', @ischar);
addParameter(p, 'groupCol',     'expcond', @ischar);
addParameter(p, 'verbose',      true,  @islogical);
parse(p, Tin, varargin{:});

thresh     = p.Results.thresh;
minN_rec   = p.Results.minN_rec;
minN_sub   = p.Results.minN_sub;
minN_group = p.Results.minN_group;
recCol     = p.Results.recCol;
subCol     = p.Results.subCol;
groupCol   = p.Results.groupCol;
logVars    = p.Results.logVars;
verbose    = p.Results.verbose;

%% ---- Identify measure columns ------------------------------------------

allVarNames = Tin.Properties.VariableNames;

if isempty(p.Results.measureCols)
    % default: all numeric columns
    isNum = varfun(@isnumeric, Tin, 'OutputFormat', 'uniform');
    measureNames = allVarNames(isNum);
else
    mc = p.Results.measureCols;
    if isnumeric(mc)
        measureNames = allVarNames(mc);
    else
        measureNames = mc;
    end
end

% exclude grouping columns from measure list even if numeric
groupingCols = {recCol, subCol, groupCol};
measureNames = measureNames(~ismember(measureNames, groupingCols));
nMeasures    = numel(measureNames);

%% ---- Standardise grouping columns to cellstr ---------------------------

Tout = Tin;   % work on a copy

for col = {recCol, subCol, groupCol}
    c = col{1};
    if ~ismember(c, allVarNames); continue; end
    if iscategorical(Tout.(c)) || ischar(Tout.(c))
        Tout.(c) = cellstr(Tout.(c));
    elseif isnumeric(Tout.(c))
        Tout.(c) = arrayfun(@num2str, Tout.(c), 'UniformOutput', false);
    end
end

n = height(Tout);

% unique grouping values
uniqueRecs   = unique(Tout.(recCol));
uniqueSubs   = unique(Tout.(subCol));
uniqueGroups = unique(Tout.(groupCol));

%% ---- Main loop over measures -------------------------------------------

report.params.thresh   = thresh;
report.params.minN_rec = minN_rec;
report.params.minN_sub = minN_sub;
report.params.minN_group = minN_group;
report.params.logVars  = logVars;
report.totalFlagged    = 0;

if verbose
    fprintf('\n==============================================\n');
    fprintf('  flagOutliersPerMeasure\n');
    fprintf('  %d measures, threshold = %.1f\n', nMeasures, thresh);
    fprintf('==============================================\n');
    fprintf('  %-35s %8s\n', 'Measure', 'Flagged');
    fprintf('  %s\n', repmat('-',1,45));
end

for mi = 1:nMeasures

    varName = measureNames{mi};

    % skip if not numeric in the output table
    if ~isnumeric(Tout.(varName))
        if verbose
            fprintf('  %-35s   SKIP (non-numeric)\n', varName);
        end
        continue
    end

    xAll     = Tout.(varName);          % full column, may contain NaN
    isLog    = ismember(varName, logVars);
    flagsVar = false(n, 1);             % will mark outliers for this measure
    levelUsed = repmat({'none'}, n, 1); % track which level flagged each row

    %% --- Per-recording pass ---
    for r = 1:numel(uniqueRecs)
        recID  = uniqueRecs(r);
        idxRec = strcmp(Tout.(recCol), recID);

        % valid = rows in this recording with non-NaN values
        valid = idxRec & ~isnan(xAll);
        x     = xAll(valid);
        if isempty(x); continue; end

        % determine reference distribution and level
        [refVals, refMask, level] = getReferenceDistribution( ...
            x, valid, xAll, Tout, subCol, groupCol, idxRec, ...
            minN_rec, minN_sub, minN_group, isLog);

        if isempty(refVals); continue; end

        % compute robust z for ALL valid rows in the reference mask
        [med, scaledMAD] = robustStats(refVals);
        if scaledMAD == 0; continue; end  % all identical — skip

        % apply to the recording rows only (not the full reference)
        idxRecValid = find(idxRec & ~isnan(xAll));
        for k = 1:numel(idxRecValid)
            row = idxRecValid(k);
            val = xAll(row);
            if isLog && val > 0
                z = abs((log10(val) - med) / scaledMAD);
            elseif isLog && val <= 0
                % can't log-transform — use raw scale with global stats
                [gMed, gMAD] = robustStats(xAll(~isnan(xAll)));
                if gMAD == 0; continue; end
                z = abs((val - gMed) / gMAD);
            else
                z = abs((val - med) / scaledMAD);
            end

            if z > thresh
                flagsVar(row)   = true;
                levelUsed{row}  = level;
            end
        end
    end % recordings

    % set outlier values to NaN
    Tout.(varName)(flagsVar) = NaN;

    nFlagged = sum(flagsVar);
    report.totalFlagged = report.totalFlagged + nFlagged;
    report.perMeasure(mi).name       = varName;
    report.perMeasure(mi).nFlagged   = nFlagged;
    report.perMeasure(mi).flaggedRows = find(flagsVar);
    report.perMeasure(mi).levelUsed  = levelUsed(flagsVar);

    if verbose
        fprintf('  %-35s %8d\n', varName, nFlagged);
    end

end % measures

if verbose
    fprintf('  %s\n', repmat('-',1,45));
    fprintf('  %-35s %8d\n', 'TOTAL (row,measure) pairs', report.totalFlagged);
    fprintf('  Rows in table: %d (unchanged)\n', height(Tout));
    fprintf('==============================================\n\n');
end

end % main function


%% =========================================================================
%  LOCAL FUNCTIONS
%% =========================================================================

function [refVals, refMask, level] = getReferenceDistribution( ...
    xRec, validRec, xAll, T, subCol, groupCol, idxRec, ...
    minN_rec, minN_sub, minN_group, isLog)
% Returns the reference distribution values, the mask of rows used,
% and a string label for which level was selected.

    refVals = []; refMask = validRec; level = 'none';

    % filter for log if needed
    if isLog
        useRec = xRec(xRec > 0);
    else
        useRec = xRec;
    end

    if numel(useRec) >= minN_rec
        if isLog; refVals = log10(useRec); else; refVals = useRec; end
        refMask = validRec;
        level   = 'recording';
        return
    end

    % subject fallback
    subjID  = unique(T.(subCol)(idxRec));
    idxSub  = ismember(T.(subCol), subjID);
    validSub = idxSub & ~isnan(xAll);
    xSub    = xAll(validSub);
    if isLog; useSub = xSub(xSub > 0); else; useSub = xSub; end

    if numel(useSub) >= minN_sub
        if isLog; refVals = log10(useSub); else; refVals = useSub; end
        refMask = validSub;
        level   = 'subject';
        return
    end

    % group fallback
    grpID   = unique(T.(groupCol)(idxRec));
    idxGrp  = ismember(T.(groupCol), grpID);
    validGrp = idxGrp & ~isnan(xAll);
    xGrp    = xAll(validGrp);
    if isLog; useGrp = xGrp(xGrp > 0); else; useGrp = xGrp; end

    if numel(useGrp) >= minN_group
        if isLog; refVals = log10(useGrp); else; refVals = useGrp; end
        refMask = validGrp;
        level   = 'group';
        return
    end

    % global fallback
    validGlob = ~isnan(xAll);
    xGlob     = xAll(validGlob);
    if isLog; useGlob = xGlob(xGlob > 0); else; useGlob = xGlob; end
    if isempty(useGlob); return; end

    if isLog; refVals = log10(useGlob); else; refVals = useGlob; end
    refMask = validGlob;
    level   = 'global';

end % getReferenceDistribution


function [med, scaledMAD] = robustStats(x)
% Compute median and scaled MAD. Falls back to std if MAD == 0.
    med = median(x);
    madv = median(abs(x - med));
    scaledMAD = 1.4826 * madv;
    if scaledMAD == 0
        s = std(x);
        if s > 0
            scaledMAD = s;
        else
            scaledMAD = 0;   % caller will skip
        end
    end
end