function vessel = getStimTrig(vessel, dN)
    % Stimulus-triggered (onset-relative) binning of the raw images, area and
    % velocity, leveraging indexTs2Trial (as in getAreaDiamVelFlowFaaProxyFromTs)
    % so we can work with stimulus-triggered data early in the pipeline. D, Q and
    % Faa are NOT computed here.
    %
    % Each onset-relative time bin keeps the pooled data points (across runs x
    % trial-columns) -- not just a summary -- so per-point analyses (incl. the later
    % D/Q/Faa) remain possible; the per-bin mean & sem are stored alongside.
    %
    % Usage:
    %   vessel = getStimTrig(vessel)        % dN=0
    %   vessel = getStimTrig(vessel, dN)
    % The design/timing (dsgn, tr, tsStartTime) is read from the vessel fields, so
    % no rCond argument is needed.
    %
    % dN : window spec (default 0). 0 -> one bin per onset-relative frame (a clean
    %      partition: each ts frame lands in exactly one bin). A scalar dN>0 gives a
    %      sliding (overlapping) window of half-width dN; [n1 n2] a single fixed
    %      window. See indexTs2Trial.
    %
    % Requires: raw images vessel.im.ts.im (per-run cell of X x Y x Z x time) and the
    % area/velocity proxies vessel.im.tsArea / tsVel (run getAreaVel with 'ts' first).
    %
    % Output: vessel(v).stimTrig with
    %   .dN              : window spec used
    %   .tStart          : 1 x nWin bin window start (s) = iStart*dt
    %   .tEnd            : 1 x nWin bin window end   (s) = (iEnd+1)*dt
    %   .t               : 1 x nWin bin window centre (s) = mean([tStart tEnd])
    %   .tStim           : 1 x nWin cell; stimulus-onset-relative time of each point
    %   .tRun            : 1 x nWin cell; run-onset-relative time of each point
    %   .trialIdx        : 1 x nWin cell; trial index within run, per point
    %   .runIdx          : 1 x nWin cell; run index, per point
    %   .raw  / .area / .vel : per-field struct with
    %       .data : 1 x nWin cell of pooled points (area/vel: 1 x nObs; raw: nVox x nObs)
    %       .mean : per-bin mean (area/vel: 1 x nWin; raw: nVox x nWin)
    %       .sem  : per-bin sem  (same shape as .mean)
    %       .n    : per-bin count of non-NaN data points (same shape as .mean)
    %   .area / .vel additionally carry the binned baseline (from .vecBase):
    %       .dataBase : 1 x nWin cell; each point's run baseline (same shape as .data),
    %                   so dA/A, dV/V can be formed later as (data-dataBase)./dataBase
    % Point ordering is shared across fields and the coordinate cells (column-major
    % over [run x trial-column]).

    if nargin<2 || isempty(dN); dN = 0; end
    % design/timing taken from the vessel fields (dsgn, tr, tsStartTime)
    rCond = rCondFromVessel(vessel);

    % build the onset-relative window index (as in getAreaDiamVelFlowFaaProxyFromTs)
    vessel = indexTs2Trial(vessel, rCond, dN);

    for v = 1:length(vessel)
        align = vessel(v).trial.align;
        dt    = align.dt;
        res   = vessel(v).trial.res(end);          % the window set just built for this dN
        nWin  = numel(res.winCols);

        % per-run data: area/vel as [nRun x nT]; raw as [nVox x nT] per run
        Marea = cat(1, vessel(v).im.tsArea.vec{:});
        Mvel  = cat(1, vessel(v).im.tsVel.vec{:});
        baseA = cellfun(@(b) b(1), vessel(v).im.tsArea.vecBase);  % 1 x nRun (one baseline per run)
        baseV = cellfun(@(b) b(1), vessel(v).im.tsVel.vecBase);
        rawC  = vessel(v).im.ts.im;
        nRun  = numel(rawC);
        szIm  = size(rawC{1}); nVox = prod(szIm(1:3));
        rawV  = cellfun(@(im) reshape(im,nVox,szIm(4)), rawC, 'UniformOutput',false);  % [nVox x nT] per run

        hw = isscalar(res.dN)*dN;                  % sliding half-width (0 for a fixed [n1 n2] window)
        % per-bin time window: a point indexed by its start time i spans [i*dt,(i+1)*dt],
        % so a bin over indices [iStart,iEnd] covers [iStart*dt, (iEnd+1)*dt]. res.tStart =
        % iStart*dt and res.tEnd = iEnd*dt, hence tEnd = res.tEnd + dt and t = window center.
        st = struct(); st.dN = res.dN;
        st.tStart = res.tStart;                    % iStart*dt
        st.tEnd   = res.tEnd + dt;                 % (iEnd+1)*dt
        st.t      = (st.tStart + st.tEnd)/2;       % window center
        [st.tStim, st.tRun, st.trialIdx, st.runIdx] = deal(cell(1,nWin));
        rawData = cell(1,nWin); areaData = cell(1,nWin); velData = cell(1,nWin);
        areaBaseData = cell(1,nWin); velBaseData = cell(1,nWin);
        rawMean = nan(nVox,nWin); rawSem = nan(nVox,nWin); rawN = nan(nVox,nWin);
        [areaMean,areaSem,areaN,velMean,velSem,velN] = deal(nan(1,nWin));

        for i = 1:nWin
            % the (trial,column) samples feeding this bin, same selection as
            % res.winCols{i}, but recovered so each point keeps its coordinates
            delay = round(res.t(i)/dt);
            mask  = ismember(align.idxTrialOnset, delay-hw : delay+hw);
            cols  = align.idxRunOnset(mask).';     % 1 x m within-run column indices
            trl   = align.trialIdx(mask).';        % 1 x m trial index
            ts    = align.idxTrialOnset(mask).';   % 1 x m onset-relative index
            m     = numel(cols);

            % coordinates expanded to [nRun x m] then flattened column-major -> 1 x (nRun*m)
            runG  = repmat((1:nRun).', 1, m);
            st.runIdx{i}   = runG(:).';
            st.trialIdx{i} = reshape(repmat(trl,       nRun,1),1,[]);
            st.tRun{i}     = reshape(repmat((cols-1)*dt,nRun,1),1,[]);
            st.tStim{i}    = reshape(repmat(ts*dt,     nRun,1),1,[]);

            % values (same column-major ordering as the coordinates)
            aBlk = Marea(:,cols); areaData{i} = aBlk(:).';
            vBlk = Mvel(:,cols);  velData{i}  = vBlk(:).';
            % per-point baseline (vecBase is one value per run -> constant across cols)
            areaBaseData{i} = reshape(repmat(baseA(:),1,m),1,[]);
            velBaseData{i}  = reshape(repmat(baseV(:),1,m),1,[]);
            R = nan(nVox,nRun,m);
            for r = 1:nRun; R(:,r,:) = rawV{r}(:,cols); end
            rawData{i} = reshape(R,nVox,nRun*m);

            % per-bin summaries (n = number of non-NaN data points contributing)
            areaMean(i) = mean(areaData{i},'omitnan'); [areaSem(i),areaN(i)] = sem1(areaData{i});
            velMean(i)  = mean(velData{i}, 'omitnan'); [velSem(i), velN(i)]  = sem1(velData{i});
            rawN(:,i)   = sum(~isnan(rawData{i}),2);
            rawMean(:,i)= mean(rawData{i},2,'omitnan'); rawSem(:,i)= std(rawData{i},0,2,'omitnan')./sqrt(rawN(:,i));
        end
        st.raw  = struct('data',{rawData}, 'mean',rawMean, 'sem',rawSem, 'n',rawN);
        st.area = struct('data',{areaData},'mean',areaMean,'sem',areaSem,'n',areaN,'dataBase',{areaBaseData});
        st.vel  = struct('data',{velData}, 'mean',velMean, 'sem',velSem, 'n',velN, 'dataBase',{velBaseData});
        vessel(v).stimTrig = st;
    end
end

function [s,n] = sem1(x)
    ok = ~isnan(x);
    n  = nnz(ok);
    s  = std(x(ok))./sqrt(n);
end

function rCond = rCondFromVessel(vessel)
    % Build an rCond-equivalent (.dsgn, .tr, .tsStartTime) from fields carried on
    % the vessel (taken from vessel(1), as one shared rCond for the array). Errors
    % if any is missing/empty. Mirrors getAreaDiamVelFlowFaaProxyFromTs.
    need    = {'dsgn','tr','tsStartTime'};
    missing = need(~isfield(vessel,need));
    for k = 1:numel(need)
        if isfield(vessel,need{k}) && isempty(vessel(1).(need{k}))
            missing{end+1} = need{k}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('getStimTrig:noRCond', ...
            ['rCond not provided and vessel is missing the fallback design/timing ' ...
             'field(s): %s. Either pass rCond (.dsgn/.tr/.tsStartTime) or set those ' ...
             'fields on vessel.'], strjoin(unique(missing),', '));
    end
    rCond             = struct();
    rCond.dsgn        = vessel(1).dsgn;
    rCond.tr          = vessel(1).tr;
    rCond.tsStartTime = vessel(1).tsStartTime;
end
