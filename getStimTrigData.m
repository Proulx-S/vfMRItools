function vessel = getStimTrigData(vessel, out, winSpec)
    % Extract stimulus-triggered data into a uniform "windows" struct by RECOMBINING
    % the finest bins from getStimTrigIdx ("binning the bins") -- no call to
    % indexTs2Trial. The result is APPENDED to vessel(v).stimTrig (a cell array), so
    % repeated calls with different windowing/data instructions accumulate rather than
    % clobber: vessel.stimTrig{1}, {2}, ... One routine, two modes, IDENTICAL output
    % structure (so windows from different calls are interchangeable):
    %   (A) timecourse : winSpec = a scalar window HALF-WIDTH (frames per side) -> overlapping,
    %                    dt-spaced windows, each 2*winSpec+1 frames wide (NOT a window count)
    %   (B) periods    : winSpec = N-by-2 [n1 n2] window bounds (n2=Inf -> last frame), or
    %                    struct('win',[N x 2],'label',{cellstr}) -> N arbitrary windows
    %   (C) all points : winSpec = inf -> ONE window pooling every frame; tStart/tEnd/t
    %                    are NaN (the window span is meaningless when all points are pooled)
    %
    %   vessel = getStimTrigData(vessel, out, winSpec)   % appends vessel.stimTrig{end+1}
    %   vessel = getStimTrigData(vessel)                 % defaults: out={'area','vel'}, winSpec=0
    %
    % winSpec : window spec (default 0 -> finest, single-frame timecourse bins); see (A)/(B)/(C).
    % out : cellstr selecting which DATA outputs to emit (default {'area','vel'}):
    %        - data signals 'raw'|'area'|'vel' -> st.<sig> = struct(.data {1xnWin cell of
    %          pooled points; area/vel 1 x nObs, raw nVox x nObs}, .mean/.sem/.n per
    %          window, and .dataBase per point per-run baseline for area/vel)
    %        - coords 'tStim'|'tRun'|'trialIdx'|'runIdx' -> st.<coord> (1xnWin cell)
    %      Window metadata is ALWAYS emitted (uniform across modes): st.winSpec, st.label,
    %      st.t, st.tStart, st.tEnd (s), st.n (1xnWin RETAINED point count, after dropping
    %      censored/NaN points; <= nRun*nCols).
    %
    % Censoring: points that are NaN in the source proxies (set upstream by applyPreproc's
    % frame censoring) are DROPPED here when the bins are formed -- so st.<sig>.data carry
    % only kept points and downstream consumers (getDiam, getDownStreamResistance, plots)
    % see no NaN and need no change. Censoring blanks whole frames, so area/vel/raw share
    % the same dropped points (one shared mask), keeping the pooled signals point-aligned.
    %
    % Requires getStimTrigIdx (vessel.stimTrigIdx) and getAreaVel ('ts'); 'raw' also
    % needs vessel.im.ts.im. Window times are taken from the constituent finest bins
    % (each frame covers [t, t+dtTs]); for winSpec they span (2*winSpec+1)*dtTs.

    if nargin<2 || isempty(out);     out     = {'area','vel'}; end
    if ischar(out); out = {out}; end
    if nargin<3 || isempty(winSpec); winSpec = 0;              end

    label0 = {};
    if isstruct(winSpec)
        win = winSpec.win; if isfield(winSpec,'label'); label0 = winSpec.label; end
        scalarSpec = false;
    elseif isscalar(winSpec)
        scalarSpec = true;
    else
        win = winSpec; scalarSpec = false;                 % N-by-2 period list
    end

    sig = intersect({'raw','area','vel'}, out);        % requested data signals
    crd = intersect({'tStim','tRun','trialIdx','runIdx'}, out);   % requested coords

    for v = 1:numel(vessel)
        assert(isfield(vessel(v),'stimTrigIdx') && ~isempty(vessel(v).stimTrigIdx), ...
            'getStimTrigData:noIdx','vessel(%d) has no stimTrigIdx (run getStimTrigIdx first).', v);
        ix = vessel(v).stimTrigIdx;
        dt = ix.dtTs;

        % --- target windows as finest-bin index ranges [lo hi] (uniform over modes) ---
        if scalarSpec
            if isinf(winSpec)
                lo = -inf; hi = inf;                              % winSpec=inf -> one window pooling all frames
            else
                d  = (min(ix.idx)+winSpec) : (max(ix.idx)-winSpec);  % valid centre delays
                lo = d - winSpec;  hi = d + winSpec;                 % each window = 2*winSpec+1 frames
            end
            label = repmat({''},1,numel(lo));
        else
            lo = win(:,1).';  hi = win(:,2).';  hi(isinf(hi)) = max(ix.idx);
            label = label0; if isempty(label); label = repmat({''},1,numel(lo)); end
        end
        nWin = numel(lo);

        % --- recombine finest bins into each window (union cols/trl/tt; span times) ---
        [winCols,winTrl,winTT] = deal(cell(1,nWin));
        [tStart,tEnd,t] = deal(nan(1,nWin));
        spanless = scalarSpec && isinf(winSpec);   % all-points pool -> window span is meaningless
        for i = 1:nWin
            b = find(ix.idx>=lo(i) & ix.idx<=hi(i));
            winCols{i} = [ix.winCols{b}];
            winTrl{i}  = [ix.winTrl{b}];
            winTT{i}   = [ix.winTT{b}];
            if spanless
                tStart(i) = NaN; tEnd(i) = NaN; t(i) = NaN;
            else
                tStart(i) = min(ix.tStart(b));
                tEnd(i)   = max(ix.tEnd(b));
                t(i)      = (tStart(i)+tEnd(i))/2;
            end
        end

        % --- source matrices for the requested signals only ---
        % Baseline proxy as [nRun x nB]: nB = nT (per-frame baseline, e.g. undetrended drift)
        % or nB = 1 (a single per-run value, e.g. the detrended global mean). baseBlock()
        % below indexes per-frame when nB>1 and broadcasts when nB==1, so dataBase pairs with
        % data correctly either way.
        if any(strcmp('area',sig)); Marea=cat(1,vessel(v).im.tsArea.vec{:}); MbaseA=cat(1,vessel(v).im.tsArea.vecBase{:}); end
        if any(strcmp('vel', sig)); Mvel =cat(1,vessel(v).im.tsVel.vec{:});  MbaseV=cat(1,vessel(v).im.tsVel.vecBase{:});  end
        if any(strcmp('raw', sig)); rawC=vessel(v).im.ts.im; szIm=size(rawC{1}); nVox=prod(szIm(1:3));
            rawV=cellfun(@(im)reshape(im,nVox,szIm(4)),rawC,'UniformOutput',false); end
        nRun = numel(vessel(v).im.tsArea.vec);

        % --- assemble (window metadata always present -> uniform structure) ---
        s = struct(); s.winSpec = winSpec; s.label = label;
        s.t = t; s.tStart = tStart; s.tEnd = tEnd; s.n = nan(1,nWin);
        for c = 1:numel(crd); s.(crd{c}) = cell(1,nWin); end
        for g = 1:numel(sig)
            s.(sig{g}).data = cell(1,nWin);
            if ~strcmp(sig{g},'raw'); s.(sig{g}).dataBase = cell(1,nWin); end
        end

        haveArea = any(strcmp('area',sig));
        haveVel  = any(strcmp('vel', sig));
        haveRaw  = any(strcmp('raw', sig));
        for i = 1:nWin
            cols = winCols{i}(:).'; m = numel(cols);
            % per-point blocks (column-major over [run x window-column])
            if haveArea; aBlk=Marea(:,cols); aDat=aBlk(:).'; aBase=baseBlock(MbaseA,cols); end
            if haveVel;  vBlk=Mvel(:,cols);  vDat=vBlk(:).'; vBase=baseBlock(MbaseV,cols); end
            if haveRaw;  R3=nan(nVox,nRun,m); for r=1:nRun; R3(:,r,:)=rawV{r}(:,cols); end; rDat=reshape(R3,nVox,nRun*m); end
            % coords (same column-major layout)
            cRun   = reshape(repmat((1:nRun).',1,m),1,[]);
            cTrial = reshape(repmat(winTrl{i}(:).',nRun,1),1,[]);
            cTRun  = reshape(repmat((cols-1)*dt,nRun,1),1,[]);
            cTStim = reshape(repmat(winTT{i}(:).',nRun,1),1,[]);
            % validity mask: drop censored/NaN points (NaN set upstream by applyPreproc) so the
            % pooled .data carry only kept points and downstream fits see no NaN. Censoring blanks
            % whole frames, so the NaN pattern is shared across area/vel/raw -> one shared mask.
            keep = true(1,nRun*m);
            if haveArea; keep = keep & isfinite(aDat); end
            if haveVel;  keep = keep & isfinite(vDat); end
            if haveRaw && ~haveArea && ~haveVel; keep = keep & ~all(isnan(rDat),1); end
            s.n(i) = nnz(keep);                       % retained (post-censoring) point count
            % store compacted (kept points only)
            if haveArea; s.area.data{i}=aDat(keep); s.area.dataBase{i}=aBase(keep); end
            if haveVel;  s.vel.data{i} =vDat(keep); s.vel.dataBase{i} =vBase(keep); end
            if haveRaw;  s.raw.data{i} =rDat(:,keep); end
            if any(strcmp('runIdx',  crd)); s.runIdx{i}   = cRun(keep);   end
            if any(strcmp('trialIdx',crd)); s.trialIdx{i} = cTrial(keep); end
            if any(strcmp('tRun',    crd)); s.tRun{i}     = cTRun(keep);  end
            if any(strcmp('tStim',   crd)); s.tStim{i}    = cTStim(keep); end
        end

        % per-window aggregates per requested signal
        for g = 1:numel(sig)
            nm = sig{g}; dd = s.(nm).data;
            if strcmp(nm,'raw')
                nv=size(dd{1},1); [mu,se,nn]=deal(nan(nv,nWin));
                for i=1:nWin; nn(:,i)=sum(~isnan(dd{i}),2); mu(:,i)=mean(dd{i},2,'omitnan'); se(:,i)=std(dd{i},0,2,'omitnan')./sqrt(nn(:,i)); end
            else
                [mu,se,nn]=deal(nan(1,nWin));
                for i=1:nWin; mu(i)=mean(dd{i},'omitnan'); [se(i),nn(i)]=sem1(dd{i}); end
            end
            s.(nm).mean=mu; s.(nm).sem=se; s.(nm).n=nn;
        end

        % append (don't clobber): accumulate this window set in the vessel.stimTrig cell
        if isfield(vessel(v),'stimTrig') && iscell(vessel(v).stimTrig) && ~isempty(vessel(v).stimTrig)
            vessel(v).stimTrig{end+1} = s;
        else
            vessel(v).stimTrig = {s};
        end
    end
end

function [s,n] = sem1(x)
    ok = ~isnan(x);
    n  = nnz(ok);
    s  = std(x(ok))./sqrt(n);
end

function blk = baseBlock(Mbase, cols)
    % Per-point baseline block (column-major over [run x window-column]) for one window.
    % Mbase is [nRun x nB]: nB>1 -> per-frame baseline, indexed at the window columns;
    % nB==1 -> a single per-run baseline value (e.g. the detrended global mean), broadcast
    % across the window columns. Returns a 1 x (nRun*numel(cols)) row to match the data.
    if size(Mbase,2)==1
        blk = repmat(Mbase, 1, numel(cols));
    else
        blk = Mbase(:,cols);
    end
    blk = blk(:).';
end
