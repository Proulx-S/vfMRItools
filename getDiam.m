function vessel = getDiam(vessel)
    % Diameter of the stimulus-triggered data, per time bin: applies computeD
    % (D = 2*sqrt(A/pi)) to each pooled area data point (and its baseline) in
    % vessel.stimTrig.area, and stores the result as vessel.stimTrig.diam, mirroring
    % the stimTrig.area struct (.data/.mean/.sem/.n/.dataBase).
    %
    % Requires getStimTrig first (vessel.stimTrig.area with .data and .dataBase).

    for v = 1:length(vessel)
        a = vessel(v).stimTrig.area;
        assert(isfield(a,'data') && isfield(a,'dataBase'), 'getDiam:noArea', ...
            'vessel(%d).stimTrig.area is missing .data/.dataBase (run getStimTrig first).', v);
        nWin = numel(a.data);

        d = struct();
        d.data     = cellfun(@computeD, a.data,     'UniformOutput', false);  % per-point diameter
        d.dataBase = cellfun(@computeD, a.dataBase, 'UniformOutput', false);  % per-point baseline diameter
        [d.mean, d.sem, d.n] = deal(nan(1,nWin));
        for i = 1:nWin
            d.mean(i) = mean(d.data{i},'omitnan');
            [d.sem(i), d.n(i)] = sem1(d.data{i});
        end
        vessel(v).stimTrig.diam = d;
    end
end

function [s,n] = sem1(x)
    ok = ~isnan(x);
    n  = nnz(ok);
    s  = std(x(ok))./sqrt(n);
end
