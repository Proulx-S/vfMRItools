function vessel = indexTs2Trial(vessel,rCond,winSpec)
    % Build the stimulus-onset-aligned trial index for a vessel's raw timeseries,
    % and the pooled time-column index sets ("windows") relative to onset. This is
    % the indexing core of getFaa.m, extracted as a reusable function so that any
    % per-(run,time) timeseries (area/velocity proxies, dX/X, dQ/Q, ...) can be
    % pooled into onset-relative windows the same way faa is -- giving every
    % windowed quantity a time axis that matches the faa timecourse.
    %
    % (Indexing logic is duplicated from getFaa.m by design; getFaa.m is left
    % untouched. Keep the two getAlign/getIdx copies in sync if either changes.)
    %
    % Usage (mirrors getFaa):
    %   vessel = indexTs2Trial(vessel, rCond{S}.(acq).(task))   % build the grid, winSpec=3
    %   vessel = indexTs2Trial(vessel, rCond{S}.(acq).(task), winSpec)
    %   vessel = indexTs2Trial(vessel, [], winSpec)                  % reuse the grid, add a window set
    %
    % Inputs:
    %   vessel : vessel struct array. Must carry the raw-timeseries area proxy
    %            (vessel(v).im.tsArea.vec, a per-run cell of [vox x time]) so the
    %            number of time points per run is known. On reuse calls (rCond=[])
    %            the stored vessel(v).trial.align is used instead.
    %   rCond  : the task struct/runCond (.dsgn onsets, .tr timeseries TR, and
    %            .tsStartTime = (nFrameOrig-nFrame)*tr, the time of the first
    %            dummy-removed ts frame rel. to the full-ts 0s start) on the first
    %            call; [] to reuse the stored grid. tsStartTime is REQUIRED (errors
    %            if missing/empty): onsetList is in full-ts time, so the ts grid is
    %            built starting at tsStartTime -> outputs are true-onset-relative.
    %   winSpec     : window spec (default 3). winSpec is the window HALF-WIDTH, NOT a count of
    %            windows.
    %            scalar  -> sliding-window timecourse: each window spans 2*winSpec+1
    %                       onset-relative frames (half-width winSpec on each side of a
    %                       centre delay); one window is produced per centre delay, and
    %                       consecutive windows are one frame (dt) apart, so they OVERLAP
    %                       for winSpec>0. winSpec=0 gives minimum-width (single-frame),
    %                       NON-overlapping windows. The forward edge never reaches the
    %                       next onset, so a CONSTANT inter-onset interval is assumed
    %                       (errors otherwise). [] -> build the grid (align) only.
    %            [n1 n2] -> a single fixed post-stim window from frame n1 to n2 (rel.
    %                       onset; n2 = inf -> up to the last frame before the next onset).
    %
    % Output: vessel, with vessel(v).trial:
    %   .dt    : timeseries sampling period (s)
    %   .align : reusable grid (dt, tt, idxTrialOnset, idxRunOnset, trialIdx)
    %   .res   : indexed struct array, one element appended per call:
    %       .winSpec      : the winSpec used for this window set
    %       .winCols : 1 x nWin cell; winCols{i} are the column indices (into the
    %                  [run x time] proxy matrices) pooled for window i
    %       .winTT   : 1 x nWin cell; onset-relative times (s) of those columns
    %       .winTrl  : 1 x nWin cell; trial index each column came from
    %       % Window time SPAN in seconds, accounting for each frame covering [t, t+dt]
    %       % (a frame's timestamp is the START of its acquisition):
    %       .tStart  : window start (s)  = (first frame index)*dt
    %       .tEnd    : window end   (s)  = (last frame index + 1)*dt
    %       .t       : window centre (s) = (tStart + tEnd)/2
    %
    % Companion functions: getFaa.m (faa from the pooled dV/V vs dD/D slope over
    % these windows) and getAreaDiamVelFlowFaaProxyFromTs.m (windowed proxies).

    if nargin<3; winSpec = 3; end

    reuse = isempty(rCond);
    if ~reuse
        dsgn = rCond.dsgn;
        tr   = rCond.tr(:);
        % timeseries sampling period (s): one shared TR across runs (the runs are
        % pooled on a common time grid, so they must share a TR)
        if ~all(abs(tr-tr(1)) < 1e-5)
            error('indexTs2Trial:nonUniformTR', ...
                'pools runs on a shared time grid; per-run TRs must match (got %s s).', ...
                mat2str(unique(tr)'));
        end
        dt = mean(tr);
        % ts time-grid offset (required, no fallback): the preprocessed (dummy-removed) ts
        % starts at tsStartTime = (nFrameOrig-nFrame)*tr in the full (with-dummy) acquisition
        % timeline, while dsgn.onsetList is in that full-ts time. The grid is built starting at
        % tsStartTime so onsets land on the true frame (the alignment, not a downstream relabel).
        hasTs = (isobject(rCond) && isprop(rCond,'tsStartTime')) || (isstruct(rCond) && isfield(rCond,'tsStartTime'));
        if ~hasTs || isempty(rCond.tsStartTime)
            error('indexTs2Trial:noTsStartTime', ['rCond.tsStartTime is missing/empty. Set it to ' ...
                '(nFrameOrig-nFrame)*tr -- the time (s) of the first preprocessed (dummy-removed) ts ' ...
                'frame relative to the full-ts 0s start. dsgn.onsetList is in full-ts time, so the ' ...
                'ts grid must start at tsStartTime.']);
        end
        tsStartTime = rCond.tsStartTime;
    end

    for v = 1:length(vessel)
        if reuse
            if ~isfield(vessel(v),'trial') || ~isfield(vessel(v).trial,'align')
                error('indexTs2Trial:noAlign', ...
                    ['called with empty rCond but vessel(%d) has no stored time grid ' ...
                     '(trial.align) to reuse. Call with rCond first.'],v);
            end
            trial = vessel(v).trial;            % reuse time grid; append a new window set
        else
            align       = getAlign(vessel(v),dsgn,dt,tsStartTime); % build time grid (full-ts time) from the design
            trial       = struct();                    % fresh; this becomes res(1)
            trial.dt          = align.dt;
            trial.tsStartTime = tsStartTime;
            trial.align       = align;
        end
        vessel(v).trial = doIt(trial,winSpec);
    end

    function align = getAlign(vessel,dsgn,dt,tsStartTime)
        % run time grid in FULL-ts time (s): the dummy-removed ts starts at tsStartTime, so
        % frame j sits at tsStartTime + (j-1)*dt. dsgn.onsetList is in full-ts time, so onsets
        % now land on the true frame and tt is true-onset-relative. (Was linspace(0,...), which
        % ignored the dummy offset and placed onsets nDummy frames late.)
        nT  = size(vessel.im.tsArea.vec{1},2);
        tts = tsStartTime + linspace(0,(nT-1)*dt,nT);
        tt  = zeros(length(dsgn.onsetList),nT);
        for i = 1:length(dsgn.onsetList)
            idx = find(ismembertol(tts,dsgn.onsetList(i),0.001));
            if isempty(idx)
                error('indexTs2Trial:onsetGrid', ...
                    ['Stimulus onset %.3fs does not fall on the timeseries grid (TR=%.4fs). ' ...
                     'Check that tr is the timeseries sampling period.'],dsgn.onsetList(i),dt);
            end
            tt(i,:) = tts - tts(idx); % time relative to stimulus onset
        end
        align.dt            = dt;
        align.tt            = tt;
        align.idxTrialOnset = round(tt./dt);
        align.idxRunOnset   = repmat(1:nT,length(dsgn.onsetList),1);
        align.trialIdx      = repmat((1:length(dsgn.onsetList))',1,nT);
    end

    function trial = doIt(trial,winSpec)
        if isempty(winSpec); return; end   % build the grid (align) only; append no window set
        align = trial.align;
        [winCols,winTT,winTrl] = getIdx(align,winSpec);
        nWin = numel(winCols);
        res = struct();
        res.winSpec      = winSpec;
        res.winCols = winCols;
        res.winTT   = winTT;
        res.winTrl  = winTrl;
        res.t       = nan(1,nWin);
        res.tStart  = nan(1,nWin);
        res.tEnd    = nan(1,nWin);
        % window time span (s): each frame's timestamp is the START of its acquisition
        % and covers [t, t+dt], so the window runs from the first frame's start to the
        % last frame's start + dt; the centre is the midpoint.
        for i = 1:nWin
            res.tStart(i) = min(winTT{i}(:));
            res.tEnd(i)   = max(winTT{i}(:)) + align.dt;
            res.t(i)      = (res.tStart(i) + res.tEnd(i))/2;
        end
        if isfield(trial,'res')
            trial.res(end+1) = res;
        else
            trial.res = res;
        end
    end

    function [winCols,winTT,winTrl] = getIdx(align,winSpec)
        % For each post-stim window position, return the time-column indices (into
        % the [run x time] proxy matrices), their onset-relative times, and the trial
        % index each column came from (winTrl, for per-point provenance/coords). Index
        % sets are run-independent: the align grid is shared across runs, so
        % pooling runs is just indexing these columns across all rows.
        isi = -unique(diff(align.idxTrialOnset(:,1),[],1)); % inter-onset interval (points)
        if ~isscalar(isi)
            error('indexTs2Trial:nonUniformISI', ...
                ['the sliding-window bound assumes a constant inter-onset ' ...
                    'interval; got ISIs %s points.'],mat2str(-isi'));
        end
        if numel(winSpec)==2
            % single fixed post-stim window [winSpec(1) winSpec(2)] (time points)
            if winSpec(2)==inf; winSpec(2) = isi-1; end
            sel     = align.idxTrialOnset>=winSpec(1) & align.idxTrialOnset<=winSpec(2);
            winCols = {align.idxRunOnset(sel)};
            winTT   = {align.tt(sel)};
            winTrl  = {align.trialIdx(sel)};
        else
            % sliding window (half-width winSpec) around each post-stim delay. Cap the
            % max delay so each trial's window stays inside its own inter-onset
            % period -- the forward edge never reaches the next onset.
            idxTrialMin   = align.idxTrialOnset(1,1)+winSpec; % run-start bound
            idxTrialMax   = isi-1-winSpec;                    % next-trial bound
            idxTrial_list = idxTrialMin:idxTrialMax;
            nDelays       = numel(idxTrial_list);
            winCols = cell(1,nDelays);
            winTT   = cell(1,nDelays);
            winTrl  = cell(1,nDelays);
            for i = 1:nDelays
                idxTrial = idxTrial_list(i);
                curIdx   = ismember(align.idxTrialOnset,idxTrial-winSpec:idxTrial+winSpec);
                cols     = align.idxRunOnset(curIdx);
                % no time column may be selected twice at one delay (it would be
                % double-counted): a repeat means trial windows overlap
                if numel(unique(cols)) < numel(cols)
                    dbstack; error('indexTs2Trial:overlapTrial','overlapping trials');
                end
                winCols{i} = cols;
                winTT{i}   = align.tt(curIdx);
                winTrl{i}  = align.trialIdx(curIdx);
            end
        end
    end
end
