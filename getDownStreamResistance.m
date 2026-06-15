function st = getDownStreamResistance(st, method)
    % Downstream-resistance index (Faa) of a stimulus-triggered windows struct, per
    % window: for each window, computeFaa fits the dV/V-vs-dD/D relationship over the
    % window's pooled points (Faa = 1/2 - 1/4*slope), with a parametric 95% CI.
    % Operates uniformly on a timecourse OR a period windows struct (st may be a struct
    % array; one element per vessel).
    %
    % method : computeFaa slope-fit method (default 'aproxSlope0' -- poly1 slope through
    %          the origin; 'aproxSlope' -- free intercept). Must be a slope method: it
    %          yields one Faa + CI per window. ('exact'/'aprox' are per-point and have no
    %          single slope/CI, so they are not valid here.)
    %
    % The fractional changes for this fit use a PER-BIN baseline -- each point is
    % referenced to its own window mean:
    %   dV/V = (vel.data  - vel.mean(win) )./vel.mean(win)
    %   dD/D = (diam.data - diam.mean(win))./diam.mean(win)
    % i.e. each window's points are centred for the within-window (across run x trial)
    % slope. NB this is distinct from the dD/D & dV/V *timecourses* shown in the
    % plots, which stay on the PER-RUN baseline (.dataBase) -- see binFrac in
    % plotSingleVesselResults.
    %
    % Requires getStimTrigData ('vel') and getDiam (st.diam).
    %
    % Output: st.faa with
    %   .val      : 1 x nWin  per-window Faa
    %   .lwr/.upr : 1 x nWin  per-window 95% CI on Faa (parametric, from the fit covariance)
    %   .n        : 1 x nWin  number of points per window (= st.vel.n)
    %   .fit      : 1 x nWin  cell of the poly1 cfit objects (dV/V vs dD/D), so the exact
    %               fit line can be reused when plotting (e.g. feval(st.faa.fit{i}, x))

    if ~exist('method','var') || isempty(method); method = 'aproxSlope0'; end

    for k = 1:numel(st)
        assert(isfield(st(k),'diam'), 'getDownStreamResistance:noDiam', ...
            'st(%d).diam is missing (run getDiam first).', k);
        nWin = numel(st(k).vel.data);

        faa = struct(); [faa.val, faa.lwr, faa.upr, faa.n] = deal(nan(1,nWin));
        faa.fit = cell(1,nWin);
        faa.n = st(k).vel.n;
        for i = 1:nWin
            switch method
                case 'aproxSlope'
                    dVoV = (st(k).vel.data{i}  - st(k).vel.dataBase{i} ) ./ st(k).vel.dataBase{i};
                    dDoD = (st(k).diam.data{i} - st(k).diam.dataBase{i}) ./ st(k).diam.dataBase{i};
                case 'aproxSlope0'
                    dVoV = (st(k).vel.data{i}  - st(k).vel.mean(i) ) ./ st(k).vel.mean(i);
                    dDoD = (st(k).diam.data{i} - st(k).diam.mean(i)) ./ st(k).diam.mean(i);
                otherwise
                    dbstack; error('code that');
            end
            if any(isnan(dVoV) | isnan(dDoD)); dbstack; error('getDownStreamResistance:nan','NaN in dVoV or dDoD'); end
            % nBoot=0 -> fast parametric CI per window (bootstrap available via nBoot>0)
            [faa.val(i),faa.lwr(i),faa.upr(i),faa.fit{i}] = computeFaa(dDoD, dVoV, method, [], 0);
        end
        st(k).faa = faa;
    end
end
