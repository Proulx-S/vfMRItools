function vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel,rCond,winSpec)
    % Windowed (onset-relative) Area/Vel/Diam, fractional-change (dX/X) and
    % volumetric-flow (dQ/Q) proxies built from the RAW timeseries, on the same
    % stimulus-onset grid as faa -- so the dX/X and dQ/Q timecourses line up in
    % time with the faa timecourse. (The companion getAreaDiamVelFlowProxy 'resp'
    % proxies instead live on the deconvolution IRF grid.) It also computes faa per
    % window via getFaa.
    %
    % How: the per-(run,time) raw-ts proxies (vessel.im.ts*, from
    % getAreaDiamVelFlowProxy) are pooled into onset-relative windows by
    % indexTs2Trial. Each window's value is the mean over its pooled (run x column)
    % data points; the SEM over that same pool is stored alongside as the per-window
    % error. Averaging the already-baseline-normalized dX/X points (rather than raw
    % proxies) keeps each run on its own .vecBase baseline, exactly as faa pools the
    % per-run dD/D & dV/V points.
    %
    % Usage (mirrors getFaa / indexTs2Trial):
    %   vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel, rCond{S}.(acq).(task))      % winSpec=3
    %   vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel, rCond{S}.(acq).(task), winSpec)
    %   vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel, [], winSpec)                     % rCond from vessel
    %
    % The design/timing struct rCond (.dsgn, .tr, .tsStartTime) may instead live on
    % the vessel itself (vessel.dsgn / vessel.tr / vessel.tsStartTime). If rCond is
    % passed it is used; if it is empty/omitted these vessel fields are used as the
    % fallback; if neither is available, this errors.
    %
    % Requires the raw-ts proxies on the input (run getAreaDiamVelFlowProxy with the
    % 'ts' field first): vessel.im.tsArea / tsVel / tsDiam and
    % tsAoA / tsVoV / tsDoD / tsQoQe / tsQoQa.
    %
    % Output: vessel, with
    %   vessel(v).fromTs.t  : 1 x nWin onset-relative window time (s) = the faa grid
    %   vessel(v).fromTs.winSpec : the window spec used
    %   vessel(v).fromTs.<name> : struct with .mean and .sem (each 1 x nWin), for
    %       <name> in Area, Vel, Diam   (raw proxies, window mean)
    %                 AoA, VoV, DoD, QoQe, QoQa  (dX/X & dQ/Q, window mean)
    %   vessel(v).faa : faa per window (see getFaa)

    if nargin<3; winSpec = 3; end
    if nargin<2; rCond = []; end

    % Resolve the design/timing struct (rCond): use the passed rCond when given,
    % otherwise fall back to the fields now carried on the vessel; error if neither
    % is available.
    if isempty(rCond)
        rCond = rCondFromVessel(vessel);
    end

    % build / reuse the onset-relative window index (sliding timecourse for winSpec scalar)
    vessel = indexTs2Trial(vessel,rCond,winSpec);

    rawFlds = {'Area','Vel','Diam'};             % raw proxies (window mean)
    frcFlds = {'AoA','VoV','DoD','QoQ'};  % fractional-change / flow proxies
    allFlds = [rawFlds frcFlds];

    for v = 1:length(vessel)
        % the window set just built for this winSpec (the sliding timecourse)
        res  = vessel(v).trial.res(end);
        nWin = numel(res.winCols);

        fromTs = struct();
        fromTs.t  = res.t;
        fromTs.winSpec = res.winSpec;
        for f = 1:numel(allFlds)
            nm   = allFlds{f};
            fldS = vessel(v).im.(['ts' nm]);
            if ~isfield(fldS,'vec') || ~iscell(fldS.vec)
                error('getAreaDiamVelFlowFaaProxyFromTs:badField', ...
                    'vessel(%d): ts%s.vec missing or not a per-run cell (run getAreaDiamVelFlowProxy with ''ts'').',v,nm);
            end
            M  = cat(1,fldS.vec{:});      % [run x time]
            mu = nan(1,nWin); se = nan(1,nWin);
            for i = 1:nWin
                pool = M(:,res.winCols{i}); pool = pool(:);
                ok   = ~isnan(pool);
                mu(i) = mean(pool(ok));
                se(i) = std(pool(ok))./sqrt(nnz(ok));
            end
            fromTs.(nm).mean = mu;
            fromTs.(nm).sem  = se;
        end
        vessel(v).fromTs = fromTs;
    end

    % faa per window from the pooled dV/V vs dD/D slope (same indices)
    vessel = getFaa(vessel);
end

function rCond = rCondFromVessel(vessel)
    % Build an rCond-equivalent (.dsgn, .tr, .tsStartTime) from the fields now
    % carried on the vessel. One shared rCond is used for the whole array (as when
    % it is passed in), so it is taken from vessel(1). Errors if any field is
    % missing/empty -- there is no rCond and no fallback.
    need    = {'dsgn','tr','tsStartTime'};
    missing = need(~isfield(vessel,need));
    for k = 1:numel(need)
        if isfield(vessel,need{k}) && isempty(vessel(1).(need{k}))
            missing{end+1} = need{k}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('getAreaDiamVelFlowFaaProxyFromTs:noRCond', ...
            ['rCond not provided and vessel is missing the fallback design/timing ' ...
             'field(s): %s. Either pass rCond (the task struct with .dsgn/.tr/.tsStartTime) ' ...
             'or set those fields on vessel.'], strjoin(unique(missing),', '));
    end
    rCond             = struct();
    rCond.dsgn        = vessel(1).dsgn;
    rCond.tr          = vessel(1).tr;
    rCond.tsStartTime = vessel(1).tsStartTime;
end
