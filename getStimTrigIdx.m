function vessel = getStimTrigIdx(vessel)
    % Build the FINEST stimulus-triggered binning: non-overlapping, minimum-width
    % (single onset-relative frame) time windows. This is the canonical index from
    % which getStimTrigData recombines ("bins the bins") into any timecourse winSpec or any
    % time-period windows. It makes the ONLY call to indexTs2Trial (with winSpec=0) and does
    % NO data extraction.
    %
    %   vessel = getStimTrigIdx(vessel)
    %
    % The design/timing (dsgn, tr, tsStartTime) is read from the vessel fields.
    %
    % Output: vessel(v).stimTrigIdx (1 x nBin, one finest bin per onset-relative frame):
    %   .dtTs            : frame period (s)
    %   .idx             : onset-relative frame index of each bin
    %   .tStart/.tEnd/.t : bin window span/centre (s) -- tStart=idx*dtTs,
    %                      tEnd=(idx+1)*dtTs, t=(idx+0.5)*dtTs (a frame covers [t,t+dtTs])
    %   .winCols         : cell; column indices (into the [run x time] proxy matrices)
    %   .winTrl          : cell; trial index each column came from
    %   .winTT           : cell; onset-relative time (s) of each column

    rCond  = rCondFromVessel(vessel);
    vessel = indexTs2Trial(vessel, rCond, 0);    % finest (winSpec=0), non-overlapping bins
    for v = 1:numel(vessel)
        r  = vessel(v).trial.res(end);
        dt = vessel(v).trial.align.dt;
        ix = struct();
        ix.dtTs    = dt;
        ix.idx     = round(r.tStart./dt);        % onset-relative frame index (tStart = idx*dt)
        ix.tStart  = r.tStart;                   % = idx*dt
        ix.tEnd    = r.tEnd;                      % = (idx+1)*dt  (set by indexTs2Trial)
        ix.t       = r.t;                         % = (idx+0.5)*dt
        ix.winCols = r.winCols;
        ix.winTrl  = r.winTrl;
        ix.winTT   = r.winTT;
        vessel(v).stimTrigIdx = ix;
    end
end

function rCond = rCondFromVessel(vessel)
    % Build an rCond-equivalent (.dsgn, .tr, .tsStartTime) from fields carried on the
    % vessel (taken from vessel(1), one shared rCond for the array). Errors if any is
    % missing/empty.
    need    = {'dsgn','tr','tsStartTime'};
    missing = need(~isfield(vessel,need));
    for k = 1:numel(need)
        if isfield(vessel,need{k}) && isempty(vessel(1).(need{k}))
            missing{end+1} = need{k}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        error('getStimTrigIdx:noRCond', ...
            ['vessel is missing the design/timing field(s): %s. Set vessel.dsgn / ' ...
             'vessel.tr / vessel.tsStartTime.'], strjoin(unique(missing),', '));
    end
    rCond             = struct();
    rCond.dsgn        = vessel(1).dsgn;
    rCond.tr          = vessel(1).tr;
    rCond.tsStartTime = vessel(1).tsStartTime;
end
