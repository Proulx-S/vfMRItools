function vessel = indexTs2Trial(vessel,rCond,dN)
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
    %   vessel = indexTs2Trial(vessel, rCond{S}.(acq).(task))   % build the grid, dN=3
    %   vessel = indexTs2Trial(vessel, rCond{S}.(acq).(task), dN)
    %   vessel = indexTs2Trial(vessel, [], dN)                  % reuse the grid, add a window set
    %
    % Inputs:
    %   vessel : vessel struct array. Must carry the raw-timeseries area proxy
    %            (vessel(v).im.tsArea.vec, a per-run cell of [vox x time]) so the
    %            number of time points per run is known. On reuse calls (rCond=[])
    %            the stored vessel(v).trial.align is used instead.
    %   rCond  : the task struct (.dsgn onsets, .tr timeseries TR) on the first
    %            call; [] to reuse the stored grid.
    %   dN     : window spec (default 3), same semantics as getFaa:
    %            scalar  -> sliding-window timecourse (half-width dN points; the
    %                       forward edge never reaches the next onset, so a CONSTANT
    %                       inter-onset interval is assumed -- errors otherwise).
    %            [n1 n2] -> a single fixed post-stim window [n1 n2] (points rel. onset;
    %                       n2 = inf -> up to the last point before the next onset).
    %
    % Output: vessel, with vessel(v).trial:
    %   .dt    : timeseries sampling period (s)
    %   .align : reusable grid (dt, tt, idxTrialOnset, idxRunOnset, trialIdx)
    %   .res   : indexed struct array, one element appended per call:
    %       .dN      : the dN used for this window set
    %       .winCols : 1 x nWin cell; winCols{i} are the column indices (into the
    %                  [run x time] proxy matrices) pooled for window i
    %       .winTT   : 1 x nWin cell; onset-relative times (s) of those columns
    %       .t       : mean   onset-relative time (s) per window
    %       .tStart  : min    onset-relative time (s) per window
    %       .tEnd    : max    onset-relative time (s) per window
    %
    % Companion functions: getFaa2.m (faa from the pooled dV/V vs dD/D slope over
    % these windows) and getAreaDiamVelFlowFaaProxyFromTs.m (windowed proxies).

    if nargin<3; dN = 3; end

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
            align       = getAlign(vessel(v),dsgn,dt); % build time grid from the design
            trial       = struct();                    % fresh; this becomes res(1)
            trial.dt    = align.dt;
            trial.align = align;
        end
        vessel(v).trial = doIt(trial,dN);
    end

    function align = getAlign(vessel,dsgn,dt)
        % run time grid (s), aligned to each stimulus onset
        nT  = size(vessel.im.tsArea.vec{1},2);
        tts = linspace(0,(nT-1)*dt,nT);
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

    function trial = doIt(trial,dN)
        align = trial.align;
        [winCols,winTT] = getIdx(align,dN);
        nWin = numel(winCols);
        res = struct();
        res.dN      = dN;
        res.winCols = winCols;
        res.winTT   = winTT;
        res.t       = nan(1,nWin);
        res.tStart  = nan(1,nWin);
        res.tEnd    = nan(1,nWin);
        for i = 1:nWin
            res.t(i)      = mean(winTT{i}(:));
            res.tStart(i) = min(winTT{i}(:));
            res.tEnd(i)   = max(winTT{i}(:));
        end
        if isfield(trial,'res')
            trial.res(end+1) = res;
        else
            trial.res = res;
        end
    end

    function [winCols,winTT] = getIdx(align,dN)
        % For each post-stim window position, return the time-column indices (into
        % the [run x time] proxy matrices) and their onset-relative times. Index
        % sets are run-independent: the align grid is shared across runs, so
        % pooling runs is just indexing these columns across all rows.
        isi = -unique(diff(align.idxTrialOnset(:,1),[],1)); % inter-onset interval (points)
        if ~isscalar(isi)
            error('indexTs2Trial:nonUniformISI', ...
                ['the sliding-window bound assumes a constant inter-onset ' ...
                    'interval; got ISIs %s points.'],mat2str(-isi'));
        end
        if numel(dN)==2
            % single fixed post-stim window [dN(1) dN(2)] (time points)
            if dN(2)==inf; dN(2) = isi-1; end
            sel     = align.idxTrialOnset>=dN(1) & align.idxTrialOnset<=dN(2);
            winCols = {align.idxRunOnset(sel)};
            winTT   = {align.tt(sel)};
        else
            % sliding window (half-width dN) around each post-stim delay. Cap the
            % max delay so each trial's window stays inside its own inter-onset
            % period -- the forward edge never reaches the next onset.
            idxTrialMin   = align.idxTrialOnset(1,1)+dN; % run-start bound
            idxTrialMax   = isi-1-dN;                    % next-trial bound
            idxTrial_list = idxTrialMin:idxTrialMax;
            nDelays       = numel(idxTrial_list);
            winCols = cell(1,nDelays);
            winTT   = cell(1,nDelays);
            for i = 1:nDelays
                idxTrial = idxTrial_list(i);
                curIdx   = ismember(align.idxTrialOnset,idxTrial-dN:idxTrial+dN);
                cols     = align.idxRunOnset(curIdx);
                % no time column may be selected twice at one delay (it would be
                % double-counted): a repeat means trial windows overlap
                if numel(unique(cols)) < numel(cols)
                    dbstack; error('indexTs2Trial:overlapTrial','overlapping trials');
                end
                winCols{i} = cols;
                winTT{i}   = align.tt(curIdx);
            end
        end
    end
end
