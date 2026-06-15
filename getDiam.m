function st = getDiam(st)
    % Diameter of a stimulus-triggered windows struct, per window: applies computeD
    % (D = 2*sqrt(A/pi)) to each pooled area data point (and its per-run baseline) in
    % st.area, and adds st.diam, mirroring st.area (.data/.dataBase/.mean/.sem/.n).
    % Operates uniformly on a timecourse OR a period windows struct (st may be a
    % struct array; one element per vessel).
    %
    % Requires getStimTrigData with 'area' in its output (st.area.data and .dataBase).

    for k = 1:numel(st)
        a = st(k).area;
        assert(isfield(a,'data') && isfield(a,'dataBase'), 'getDiam:noArea', ...
            'st(%d).area is missing .data/.dataBase (request ''area'' in getStimTrigData).', k);
        nWin = numel(a.data);

        d = struct();
        d.data     = cellfun(@computeD, a.data,     'UniformOutput', false);  % per-point diameter
        d.dataBase = cellfun(@computeD, a.dataBase, 'UniformOutput', false);  % per-point baseline diameter
        [d.mean, d.sem, d.n] = deal(nan(1,nWin));
        for i = 1:nWin
            d.mean(i) = mean(d.data{i},'omitnan');
            [d.sem(i), d.n(i)] = sem1(d.data{i});
        end
        st(k).diam = d;
    end
end

function [s,n] = sem1(x)
    ok = ~isnan(x);
    n  = nnz(ok);
    s  = std(x(ok))./sqrt(n);
end
