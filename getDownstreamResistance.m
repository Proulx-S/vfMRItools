function vessel = getDownStreamResistance(vessel)
    % Downstream-resistance index (Faa) of the stimulus-triggered data, per time
    % bin: for each bin, computeFaa fits the dV/V-vs-dD/D slope through the origin
    % over the bin's pooled points (method 'aproxSlope0': Faa = 1/2 - 1/4*slope).
    %
    % The fractional changes for this fit use a PER-BIN baseline -- each point is
    % referenced to its own bin mean:
    %   dV/V = (vel.data  - vel.mean(bin) )./vel.mean(bin)
    %   dD/D = (diam.data - diam.mean(bin))./diam.mean(bin)
    % i.e. each bin's points are centred for the within-bin (across run x trial)
    % slope. NB this is distinct from the dD/D & dV/V *timecourses* shown in the
    % plots, which stay on the PER-RUN baseline (.dataBase) -- see binFrac in
    % plotSingleVesselResults.
    %
    % Requires getStimTrig (vessel.stimTrig.vel) and getDiam (vessel.stimTrig.diam).
    %
    % Output: vessel.stimTrig.faa with
    %   .val      : 1 x nWin  per-bin Faa
    %   .lwr/.upr : 1 x nWin  per-bin 95% CI on Faa (parametric, from the fit covariance)
    %   .n        : 1 x nWin  number of points per bin (= stimTrig.vel.n)

    for v = 1:length(vessel)
        st = vessel(v).stimTrig;
        assert(isfield(st,'diam'), 'getDownStreamResistance:noDiam', ...
            'vessel(%d).stimTrig.diam is missing (run getDiam first).', v);
        nWin = numel(st.vel.data);

        faa   = struct(); [faa.val, faa.lwr, faa.upr, faa.n] = deal(nan(1,nWin));
        faa.n = st.vel.n;
        for i = 1:nWin
            dVoV = (st.vel.data{i}  - st.vel.mean(i) ) ./ st.vel.mean(i);
            dDoD = (st.diam.data{i} - st.diam.mean(i)) ./ st.diam.mean(i);
            % dVoV = (st.vel.data{i}  - st.vel.dataBase{i} ) ./ st.vel.dataBase{i};
            % dDoD = (st.diam.data{i} - st.diam.dataBase{i}) ./ st.diam.dataBase{i};
            if any(isnan(dVoV) | isnan(dDoD)); dbstack; error('getDownStreamResistance:nan','NaN in dVoV or dDoD'); end
            % nBoot=0 -> fast parametric CI per bin (bootstrap available via nBoot>0)
            [faa.val(i),faa.lwr(i),faa.upr(i)] = computeFaa(dDoD, dVoV, 'aproxSlope0', [], 0);
        end
        vessel(v).stimTrig.faa = faa;
    end
end
