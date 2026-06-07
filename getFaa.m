function vessel = getFaa(vessel,rCond,dN)
    % faa index per vessel, from the area/velocity
    % proxy timeseries produced by getAreaDiamVelProxy.
    %
    % faa = 1/2 - 1/4 * (slope of dV/V vs dD/D), the poly1 fit of dV/V
    % against dD/D. dD/D is derived
    % from the area proxy via D = 2*sqrt(A/pi).
    %
    % Usage:
    %   vessel  = getFaa(vessel, rCond{S}.(acq).(task))   % builds the time grid
    %   vessel2 = getFaa(vessel, [], dN)                  % reuse the time grid
    %
    % First call: dsgn (stimulus onsets) and the timeseries TR are read from the
    % rCond task struct (.dsgn and .tr). The TR feeds the time grid used to align
    % the timeseries to stimulus onsets, so it is the timeseries sampling period
    % (acquisition TR, seconds), not the deconvolution TR. The resulting time
    % grid is stored in vessel(v).faa.align.
    %
    % Subsequent calls: pass rCond = [] and the time grid is reused from
    % vessel(v).faa.align (no rCond needed). The input vessel must already carry
    % .faa.align from a previous getFaa call.
    %
    % dN (optional, default 3): post-stim window, in time points.
    %   scalar  -> sliding-window timecourse; dN is the half-window (points on
    %              each side of each post-stim delay). The timecourse extends
    %              only as far as keeps each trial's window within its own
    %              inter-stimulus period -- the forward edge never reaches the
    %              next stimulus onset.
    %              NOTE: the sliding-window max delay is derived from a single
    %              inter-onset interval, so a CONSTANT ISI across trials is
    %              assumed. getFaa errors (getFaa:nonUniformISI) if the onsets
    %              in dsgn are not equally spaced.
    %   [n1 n2] -> a single faa pooled over the fixed post-stim window from
    %              time point n1 to n2 (relative to stimulus onset).
    %
    % Output: vessel, with vessel(v).faa:
    %   .dt     : timeseries sampling period (s)
    %   .align  : reusable time grid (dt, tt, idxTrialOnset, idxRunOnset)
    %   .all    : scalar faa using all time points pooled together (call-independent)
    %   .res    : indexed struct array, one element appended per getFaa call
    %             (the rCond call creates res(1); each reuse call appends res(end+1)):
    %       .dN     : the dN used for this result
    %       .ts     : faa value(s) -- [1 x nDelay] for the sliding-window
    %                 timecourse, scalar for a fixed [n1 n2] window
    %       .t      : mean post-stim time (s) for each .ts entry
    %       .tStart : window start time (s) for each .ts entry
    %       .tEnd   : window end time (s) for each .ts entry
    %
    % Mirrors the "faa using all time points", "faa timecourse" and
    % "faa during/post stim" sections of projects/dVdA/doIt.m.

    if nargin<3; dN = 3; end

    reuse = isempty(rCond);
    if ~reuse
        dsgn = rCond.dsgn;
        tr   = rCond.tr(:);
        % timeseries sampling period (s): one shared TR across runs (the runs
        % are pooled on a common time grid, so they must share a TR)
        if ~all(abs(tr-tr(1)) < 1e-5)
            error('getFaa:nonUniformTR', ...
                'getFaa pools runs on a shared time grid; per-run TRs must match (got %s s).', ...
                mat2str(unique(tr)'));
        end
        dt = mean(tr);
    end

    for v = 1:length(vessel)
        if reuse
            if ~isfield(vessel(v),'faa') || ~isfield(vessel(v).faa,'align')
                error('getFaa:noAlign', ...
                    ['getFaa called with empty rCond but vessel(%d) has no stored time grid ' ...
                     '(faa.align) to reuse. Run getFaa with rCond first.'],v);
            end
            faa = vessel(v).faa;                 % reuse time grid; append a new result
        else
            align      = getAlign(vessel(v),dsgn,dt); % build time grid from the design
            faa        = struct();                    % fresh; this result becomes res(1)
            faa.dt     = align.dt;
            faa.align  = align;
        end

        vessel(v).faa = doIt(vessel(v),faa,dN);
        % winCols{v} = vessel(v).faa.winCols;
    end
    % save winCols

    function align = getAlign(vessel,dsgn,dt)
        % run time grid (s), aligned to each stimulus onset
        nT  = size(vessel.im.tsArea.vec{1},2);
        tts = linspace(0,(nT-1)*dt,nT);
        tt  = zeros(length(dsgn.onsetList),nT);
        for i = 1:length(dsgn.onsetList)
            idx = find(ismembertol(tts,dsgn.onsetList(i),0.001));
            if isempty(idx)
                error('getFaa:onsetGrid', ...
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

    function faa = doIt(vessel,faa,dN)
        align = faa.align;       % reusable time grid (see getAlign)

        % per-run proxies -> [run x time], gathered across runs before fitting
        % Ats = cat(1,vessel.im.tsArea.vec{:}); Ats(Ats<0) = nan; % area (drop <0 -> imaginary D)
        Vts = cat(1,vessel.im.tsVel.vec{:});                    % velocity
        Dts = cat(1,vessel.im.tsD.vec{:});                      % diameter

        % fractional change relative to each run's mean
        dDoDts = (Dts - mean(Dts,2,'omitnan')) ./ mean(Dts,2,'omitnan');
        dVoVts = (Vts - mean(Vts,2)) ./ mean(Vts,2);

        % faa using all time points, all runs pooled (call-independent)
        if ~isfield(faa,'all') || isempty(faa.all)
            faa.all = fitFaa(dDoDts(:),dVoVts(:));
        end

        % time-column index sets for each post-stim window (run-independent)
        [winCols,winTT] = getIdx(align,dN);

        res.dN     = dN;
        nWin       = numel(winCols);
        res.ts     = nan(1,nWin);
        res.t      = nan(1,nWin);
        res.tStart = nan(1,nWin);
        res.tEnd   = nan(1,nWin);
        for i = 1:nWin
            cols = winCols{i};
            % pool runs (rows) x window columns, then slope-only fit
            res.ts(i)     = fitFaa(dDoDts(:,cols),dVoVts(:,cols));
            res.t(i)      = mean(winTT{i}(:));
            res.tStart(i) = min(winTT{i}(:));
            res.tEnd(i)   = max(winTT{i}(:));
        end

        % append this result to the indexed substructure
        if isfield(faa,'res')
            faa.res(end+1) = res;
        else
            faa.res = res;
        end
        faa.winCols = winCols;
    end

    function [winCols,winTT] = getIdx(align,dN)
        % For each post-stim window position, return the time-column indices
        % (into the [run x time] proxy matrices) and their onset-relative
        % times. Index sets are run-independent: the align grid is shared
        % across runs, so fitting pools runs by indexing these columns.
        isi = -unique(diff(align.idxTrialOnset(:,1),[],1)); % inter-onset interval (points)
        if ~isscalar(isi)
            error('getFaa:nonUniformISI', ...
                ['getFaa''s sliding-window bound assumes a constant inter-onset ' ...
                    'interval; got ISIs %s points.'],mat2str(-isi'));
        end
        if numel(dN)==2
            % single fixed post-stim window [dN(1) dN(2)] (time points)
            if dN(2)==inf; dN(2) = isi-1; end
            sel     = align.idxTrialOnset>=dN(1) & align.idxTrialOnset<=dN(2);
            winCols = {align.idxRunOnset(sel)};
            winTT   = {align.tt(sel)};
        else
            % faa timecourse (sliding window, half-width dN, around each
            % post-stim delay). Cap the max delay so each trial's window stays
            % inside its own inter-onset period -- the forward edge never
            % reaches the next onset. NOTE: a CONSTANT inter-onset interval is
            % assumed (see header).
            idxTrialMin = align.idxTrialOnset(1,1)+dN; % run-start bound
            idxTrialMax   = isi-1-dN; % next-trial bound
            idxTrial_list = idxTrialMin:idxTrialMax;
            nDelays       = numel(idxTrial_list);
            winCols = cell(1,nDelays);
            winTT   = cell(1,nDelays);
            for i = 1:nDelays
                idxTrial = idxTrial_list(i);
                curIdx   = ismember(align.idxTrialOnset,idxTrial-dN:idxTrial+dN);
                cols     = align.idxRunOnset(curIdx);
                % no time column may be selected twice at one delay (it would be
                % double-counted in the fit): a repeat means trial windows overlap
                if numel(unique(cols)) < numel(cols)
                    dbstack; error('getFaa:overlapTrial','overlapping trials');
                end
                winCols{i} = cols;
                winTT{i}   = align.tt(curIdx);
            end
        end
    end

    function faa = fitFaa(X,Y)
        ok = ~isnan(X(:)) & ~isnan(Y(:)); % drop NaNs (e.g. area<0 patched above)
        % f  = fit(X(ok),Y(ok),fittype({'x'})); % slope-only fit (intercept fixed at 0)
        % faa = 1/2 - 1/4*f.a;
        f  = fit(X(ok),Y(ok),'poly1');
        faa = 1/2 - 1/4*f.p1;
    end
end
