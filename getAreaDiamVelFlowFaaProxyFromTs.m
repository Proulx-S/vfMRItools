function vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel,rCond,dN)
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
    %   vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel, rCond{S}.(acq).(task))      % dN=3
    %   vessel = getAreaDiamVelFlowFaaProxyFromTs(vessel, rCond{S}.(acq).(task), dN)
    %
    % Requires the raw-ts proxies on the input (run getAreaDiamVelFlowProxy with the
    % 'ts' field first): vessel.im.tsArea / tsVel / tsDiam and
    % tsAoA / tsVoV / tsDoD / tsQoQe / tsQoQa.
    %
    % Output: vessel, with
    %   vessel(v).fromTs.t  : 1 x nWin onset-relative window time (s) = the faa grid
    %   vessel(v).fromTs.dN : the window spec used
    %   vessel(v).fromTs.<name> : struct with .mean and .sem (each 1 x nWin), for
    %       <name> in Area, Vel, Diam   (raw proxies, window mean)
    %                 AoA, VoV, DoD, QoQe, QoQa  (dX/X & dQ/Q, window mean)
    %   vessel(v).faa : faa per window (see getFaa)

    if nargin<3; dN = 3; end

    % build / reuse the onset-relative window index (sliding timecourse for dN scalar)
    vessel = indexTs2Trial(vessel,rCond,dN);

    rawFlds = {'Area','Vel','Diam'};             % raw proxies (window mean)
    frcFlds = {'AoA','VoV','DoD','QoQe','QoQa'};  % fractional-change / flow proxies
    allFlds = [rawFlds frcFlds];

    for v = 1:length(vessel)
        % the window set just built for this dN (the sliding timecourse)
        res  = vessel(v).trial.res(end);
        nWin = numel(res.winCols);

        fromTs = struct();
        fromTs.t  = res.t;
        fromTs.dN = res.dN;
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
