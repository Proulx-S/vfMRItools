function vessel = applyPreproc(vessel, detrendFlag, censorFlag)
    % Source-timeseries preprocessing, applied IN PLACE to vessel.im.ts before the
    % area/velocity proxies are formed (getAreaVel). Occupies the pipeline slot that
    % getTsBase used to hold and calls getTsBase itself, then optionally censors and
    % detrends (in that order):
    %
    %   censorFlag  (default 0 = none): blank out (set to NaN) the censored frames of
    %       vessel.im.ts, using one of the per-run censor masks built by getCnsr:
    %         0 | false | '' | 'none'  -> no censoring
    %         'cnsr'                    -> vessel.im.cnsr            (.vec, true = exclude)
    %         'cnsr_mainClust'          -> vessel.im.cnsr_mainClust
    %       Censored frames become NaN; getStimTrigData drops NaN points when it forms
    %       the bins, so downstream (getDiam / getDownStreamResistance / plots) needs no
    %       change. Applied BEFORE detrending (a NaN frame stays NaN through detrending).
    %
    %   detrendFlag (default false): detrend each raw point -- subtract its own full
    %       polynomial baseline (vessel.im.tsBase, from getTsBase) and add back the
    %       within-run mean baseline. Removes slow drift while preserving the run
    %       signal level. (Moved here from getAreaVel.)
    %
    % Both steps modify vessel.im.ts.im IN PLACE, so getAreaVel(...,'ts',...) consumes the
    % preprocessed timeseries. The dedicated baseline (vessel.im.basePolyRun) is left
    % untouched, so the baseline proxies (dA/A, dV/V denominators) stay drift-free and
    % uncensored.
    %
    %   vessel = applyPreproc(vessel, detrendFlag, censorFlag)

    if nargin<2 || isempty(detrendFlag); detrendFlag = false; end
    if nargin<3 || isempty(censorFlag);  censorFlag  = 0;     end
    censorField = resolveCensor(censorFlag);

    vessel = getTsBase(vessel);     % reconstruct the full polynomial-baseline timeseries (always)
    if ~isfield(vessel.im.ts,'info'); vessel.im.ts.info = {}; end

    

    for v = 1:numel(vessel)
        nRun = numel(vessel(v).im.ts.im);

        % --- censor: blank excluded frames to NaN (before detrend) -------------------
        if ~isempty(censorField)
            assert(isfield(vessel(v).im,censorField) && ~isempty(vessel(v).im.(censorField)), ...
                'applyPreproc:noCensor','vessel(%d).im.%s is missing (run getCnsr first).', v, censorField);
            cn = vessel(v).im.(censorField);
            assert(numel(cn.vec)==nRun, 'applyPreproc:censorRun', ...
                'vessel(%d).im.%s has %d runs but ts has %d.', v, censorField, numel(cn.vec), nRun);
            for r = 1:nRun
                bad = logical(cn.vec{r}(:).');                  % 1 x nFrame, true = exclude
                assert(numel(bad)==size(vessel(v).im.ts.im{r},4), 'applyPreproc:censorLen', ...
                    'vessel(%d) run %d: censor has %d frames, ts has %d.', ...
                    v, r, numel(bad), size(vessel(v).im.ts.im{r},4));
                vessel(v).im.ts.im{r}(:,:,:,bad)     = NaN;         % blank censored frames -> NaN
                vessel(v).im.tsBase.im{r}(:,:,:,bad) = NaN;         % blank censored frames -> NaN
            end
            vessel.im.ts.info{end+1} = 'censored';
        end

        % --- detrend: remove slow drift, keep the run signal level -------------------
        if detrendFlag
            assert(isfield(vessel(v).im,'tsBase') && ~isempty(vessel(v).im.tsBase), ...
                'applyPreproc:noBaseTs','detrend requires vessel.im.tsBase (from getTsBase).');
            bm = mean(cat(4,vessel(v).im.tsBase.im{:}),4,'omitnan'); 
            for r = 1:nRun
                vessel(v).im.ts.im{r}     = vessel(v).im.ts.im{r} - vessel(v).im.tsBase.im{r} + bm;
                vessel(v).im.tsBase.im{r} = bm;
            end
            vessel.im.ts.info{end+1} = 'detrended';
        end
    end
end

function fld = resolveCensor(flag)
    % Map censorFlag to the vessel.im field name to use ('' = no censoring).
    if isempty(flag) || (islogical(flag) && ~flag) || (isnumeric(flag) && isscalar(flag) && flag==0)
        fld = ''; return;
    end
    if ischar(flag) || isstring(flag)
        fld = char(flag);
        if strcmpi(fld,'none'); fld = ''; return; end
        if ~ismember(fld,{'cnsr','cnsr_mainClust'})
            error('applyPreproc:censorFlag', ...
                'censorFlag must be 0/''none'', ''cnsr'' or ''cnsr_mainClust'' (got ''%s'').', fld);
        end
        return;
    end
    error('applyPreproc:censorFlag','censorFlag must be 0/false/''none'', ''cnsr'' or ''cnsr_mainClust''.');
end
