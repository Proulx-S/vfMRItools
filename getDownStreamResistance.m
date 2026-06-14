function vessel = getDownStreamResistance(vessel)
    % Downstream-resistance index (Faa) of the stimulus-triggered data, per time
    % bin: for each bin, computeFaa fits the dV/V-vs-dD/D slope over the bin's pooled
    % points (method 'aproxSlope': Faa = 1/2 - 1/4*slope). The fractional changes are
    % formed per point from the stimTrig velocity and diameter and their baselines:
    %   dV/V = (vel.data  - vel.dataBase )./vel.dataBase
    %   dD/D = (diam.data - diam.dataBase)./diam.dataBase
    %
    % Requires getStimTrig (vessel.stimTrig.vel) and getDiam (vessel.stimTrig.diam).
    %
    % Output: vessel.stimTrig.faa with
    %   .val : 1 x nWin  per-bin Faa (NaN where < 2 valid points)
    %   .n   : 1 x nWin  number of valid (non-NaN) point pairs used per bin

    for v = 1:length(vessel)
        st = vessel(v).stimTrig;
        assert(isfield(st,'diam'), 'getDownStreamResistance:noDiam', ...
            'vessel(%d).stimTrig.diam is missing (run getDiam first).', v);
        nWin = numel(st.vel.data);

        faa = struct(); [faa.val, faa.n] = deal(nan(1,nWin));
        for i = 1:nWin
            dVoV = (st.vel.data{i}  - st.vel.dataBase{i} ) ./ st.vel.dataBase{i};
            dDoD = (st.diam.data{i} - st.diam.dataBase{i}) ./ st.diam.dataBase{i};
            ok   = ~isnan(dVoV) & ~isnan(dDoD);
            faa.n(i) = nnz(ok);
            if faa.n(i) >= 2
                faa.val(i) = computeFaa(dDoD(ok), dVoV(ok), 'aproxSlope0');
            end
        end
        vessel(v).stimTrig.faa = faa;
    end
end
