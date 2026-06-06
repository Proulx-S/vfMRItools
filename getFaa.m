function vessel = getFaa(vessel,rCond,dN)
    % Flow-area-adjustment index (faa) per vessel, from the area/velocity
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
    %              next stimulus onset (bounded by the run end and the smallest
    %              inter-onset gap).
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
    end

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
    end

    function faa = doIt(vessel,faa,dN)
        align         = faa.align;       % reusable time grid (see getAlign)
        tt            = align.tt;
        idxTrialOnset = align.idxTrialOnset;
        idxRunOnset   = align.idxRunOnset;

        % per-run proxies -> [run x time]
        % Ats = cat(1,vessel.im.tsArea.vec{:}); Ats(Ats<0) = nan; % area (drop <0 -> imaginary D)
        Vts = cat(1,vessel.im.tsVel.vec{:});                    % velocity
        Dts = cat(1,vessel.im.tsD.vec{:});                      % diameter

        % fractional change relative to each run's mean
        dDoDts = (Dts - mean(Dts,2)) ./ mean(Dts,2);
        dVoVts = (Vts - mean(Vts,2)) ./ mean(Vts,2);

        % faa using all time points (call-independent)
        faa.all = fitFaa(dDoDts(:),dVoVts(:));

        % windowed result for this dN
        res.dN = dN;
        if numel(dN)==2
            % faa within a fixed post-stim window [dN(1) dN(2)] (time points)
            sel = idxTrialOnset>=dN(1) & idxTrialOnset<=dN(2);
            idx = idxRunOnset(sel);
            T   = tt(sel);
            res.ts     = fitFaa(dDoDts(:,idx),dVoVts(:,idx));
            res.t      = mean(T(:));
            res.tStart = min(T(:));
            res.tEnd   = max(T(:));
        else
            % faa timecourse (sliding window around each post-stim delay).
            % Cap the max delay so each trial's window stays inside its own
            % period: bound nnMax by the run end AND by the inter-onset gap, so
            % the window's forward edge (nn+dN) never reaches the next onset.
            nnMax      = min(idxRunOnset(:,end) - idxRunOnset(idxTrialOnset==dN)); % run-end bound
            onsetIdx   = idxRunOnset(idxTrialOnset==0);                            % onset index per trial
            if numel(onsetIdx) > 1
                isiPts = min(diff(sort(onsetIdx(:))));   % smallest inter-onset gap (points)
                nnMax  = min(nnMax, isiPts-dN-1);        % keep nn+dN < next onset
            end
            nnMax      = max(nnMax,0);                    % guard against tiny periods
            res.ts     = nan(1,nnMax+1);
            res.t      = nan(1,nnMax+1);
            res.tStart = nan(1,nnMax+1);
            res.tEnd   = nan(1,nnMax+1);
            for nn = 0:nnMax % loop over time points after stimulus onset
                % Pool data points for the current post-stimulus delay
                idxT = idxRunOnset(idxTrialOnset==0) + nn;
                idx  = idxT;
                T    = tt(sub2ind(size(tt), 1:size(tt,1), idxT(:)'));
                for i = 1:dN
                    idxT = idxRunOnset(idxTrialOnset==-i) + nn;
                    idx  = [idx; idxT];
                    T    = [T tt(sub2ind(size(tt), 1:size(tt,1), idxT(:)'))];
                    idxT = idxRunOnset(idxTrialOnset== i) + nn;
                    idx  = [idx; idxT];
                    T    = [T tt(sub2ind(size(tt), 1:size(tt,1), idxT(:)'))];
                end
                res.ts(nn+1)     = fitFaa(dDoDts(:,idx),dVoVts(:,idx));
                res.t(nn+1)      = mean(T(:));
                res.tStart(nn+1) = min(T(:));
                res.tEnd(nn+1)   = max(T(:));
            end
        end

        % append this result to the indexed substructure
        if isfield(faa,'res')
            faa.res(end+1) = res;
        else
            faa.res = res;
        end
    end

    function faa = fitFaa(X,Y)
        % ok = ~isnan(X(:)) & ~isnan(Y(:)); % drop NaNs (e.g. area<0 patched above)
        % f  = fit(X(ok),Y(ok),fittype({'x'}));
        f  = fit(X(:),Y(:),fittype({'x'})); % slope-only fit (intercept fixed at 0)
        faa = 1/2 - 1/4*f.a;
    end
end
