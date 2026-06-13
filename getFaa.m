function vessel = getFaa(vessel)
    % faa from the dV/V vs dD/D poly1 slope, computed per onset-relative window
    % using the trial index built by indexTs2Trial (vessel.trial). The fit pools
    % the (run x window-column) dD/D & dV/V points for each window, within the
    % from-ts proxy pipeline (see getAreaDiamVelFlowFaaProxyFromTs / indexTs2Trial).
    %
    % faa = 1/2 - 1/4 * (slope of dV/V vs dD/D).
    %
    % (Formerly getFaa2; it replaced the older onset-aligning getFaa, which built its
    %  own 0-based time grid. This version inherits the tsStartTime-correct grid from
    %  indexTs2Trial, so no dummy-onset correction is needed downstream.)
    %
    % Requires, per vessel:
    %   vessel(v).trial      : from indexTs2Trial (.align, .res with .winCols/.t/...)
    %   vessel(v).im.tsDoD   : dD/D per-run timeseries (.vec, per-run cell)
    %   vessel(v).im.tsVoV   : dV/V per-run timeseries (.vec, per-run cell)
    %
    % Output: vessel, with vessel(v).faa:
    %   .dt, .align
    %   .all / .allYint / .allXint : faa & intercepts pooling all time points
    %   .res(k) : one per trial.res(k) window set, with
    %       .dN, .ts, .yint, .xint, .t, .tStart, .tEnd  (.ts is [1 x nWin])

    for v = 1:length(vessel)
        if ~isfield(vessel(v),'trial') || ~isfield(vessel(v).trial,'res')
            error('getFaa:noTrial','vessel(%d): run indexTs2Trial first.',v);
        end

        % per-run dD/D, dV/V pooled across runs -> [run x time] before fitting
        dDoDts = cat(1,vessel(v).im.tsDoD.vec{:});
        dVoVts = cat(1,vessel(v).im.tsVoV.vec{:});

        faaS = struct();
        faaS.dt    = vessel(v).trial.dt;
        faaS.align = vessel(v).trial.align;

        % faa using all time points, all runs pooled (window-independent)
        faaS.all = computeFaa(mean(dDoDts(:),[1 2]),mean(dVoVts(:),[1 2]),'exact');
        % [faaS.all,faaS.allYint,faaS.allXint] = fitFaa(dDoDts(:),dVoVts(:));

        for k = 1:numel(vessel(v).trial.res)
            tr   = vessel(v).trial.res(k);
            nWin = numel(tr.winCols);
            res = struct();
            res.dN     = tr.dN;
            res.Faa    = nan(1,nWin);
            res.yint   = nan(1,nWin);
            res.xint   = nan(1,nWin);
            res.t      = tr.t;
            res.tStart = tr.tStart;
            res.tEnd   = tr.tEnd;
            for i = 1:nWin
                cols = tr.winCols{i};
                % pool runs (rows) x window columns, then poly1 fit
                % res.Faa(i) = median(computeFaa(dDoDts(:,cols),dVoVts(:,cols),'exact'),[1 2]);
                % res.Faa(i) = computeFaa(median(dDoDts(:,cols),[1 2]),median(dVoVts(:,cols),[1 2]),'exact');
                % res.Faa(i) = computeFaa(median(dDoDts(:,cols),[1 2]),median(dVoVts(:,cols),[1 2]),'aprox');
                % res.Faa(i) = median(computeFaa(dDoDts(:,cols),dVoVts(:,cols),'aprox'),[1 2]);
                [res.Faa(i),res.yint(i),res.xint(i)] = computeFaa(dDoDts(:,cols),dVoVts(:,cols),'aproxSlope');
                % [res.Faa(i),res.yint(i),res.xint(i)] = computeFaa(dDoDts(:,cols),dVoVts(:,cols),'aproxSlope0');
                
                
                % [res.ts(i),res.yint(i),res.xint(i)] = fitFaa(dDoDts(:,cols),dVoVts(:,cols));
            end
            if isfield(faaS,'res'); faaS.res(end+1) = res; else; faaS.res = res; end
        end

        vessel(v).faa = faaS;
    end

    % function [faa,yint,xint] = fitFaa(X,Y)
    %     ok = ~isnan(X(:)) & ~isnan(Y(:)); % drop NaNs (e.g. area<0 patched upstream)
    %     f    = fit(X(ok),Y(ok),'poly1');
    %     faa  = 1/2 - 1/4*f.p1;
    %     yint = f.p2;          % dV/V at dD/D=0 (fit y-intercept)
    %     xint = -f.p2./f.p1;   % dD/D at dV/V=0 (fit x-intercept)
    % end
end
