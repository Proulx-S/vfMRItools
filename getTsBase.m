function vessel = getTsBase(vessel)
    % Reconstruct the full polynomial-baseline timeseries per run from the
    % 3dDeconvolve outputs, stored on the vessel ROI grid in the same struct
    % format as vessel.im.ts (so it can be processed like any other source field).
    %
    % For each vessel the 3dDeconvolve file triplet is derived from
    % vessel.im.resp.fName (cond-stim respAv -> cond-FULL stats):
    %   <prefix>+orig.BRIK / .HEAD : stats bucket; the sub-bricks labelled
    %        Run#rPol#p_Coef hold the per-run polynomial (baseline) coefficients
    %   <prefix>.xmat.1D           : the 3dDeconvolve design matrix; its baseline
    %        (polynomial) columns share those labels (minus '_Coef')
    % The full baseline timeseries is, per voxel, X_baseline * beta_baseline; it is
    % split per run (xmat RunStart) and cropped to the ROI exactly as vessel.im.ts
    % was. Its run-mean equals vessel.im.basePolyRun.
    %
    % Baseline coefficients / columns are selected by label (not by position), and
    % matched to each other by name, so the reconstruction is robust to bucket
    % layout, run count and polynomial order.
    %
    % AFNI is used to read the +orig bucket and the .xmat.1D (system calls; the
    % +orig dataset is referenced WITHOUT the .BRIK/.HEAD suffix). Requires the
    % global src.afni module-load prefix (set in the doIt env-setup section).
    %
    % Output: vessel(v).im.tsBase, same layout as vessel(v).im.ts
    %   (.fName, .x, .y, .mask, .im) with .im a per-run cell of [X x Y x 1 x time].

    global src %#ok<GVMIS>
    if isempty(src) || ~isfield(src,'afni') || isempty(src.afni)
        error('getTsBase:noAfni', ...
            'global src.afni (AFNI module-load prefix) is required, as set in the doIt env-setup section.');
    end

    for v = 1:length(vessel)
        vessel(v).im.tsBase = doIt(vessel(v), src.afni);
    end

    function tsBase = doIt(vessel, afni)
        % 3dDeconvolve file triplet from the resp filename
        rf = vessel.im.resp.fName; if iscell(rf); rf = rf{1}; end
        bp   = replace(replace(rf,'_cond-stim_','_cond-FULL_'),'_respAv.nii.gz','_stats');
        dset = [bp '+orig'];               % BRIK/HEAD pair, referenced without the suffix
        xmat = [bp '.xmat.1D'];
        assert(exist([dset '.BRIK'],'file')>0 && exist([dset '.HEAD'],'file')>0, ...
            'getTsBase:noStats','missing 3dDeconvolve stats %s (.BRIK/.HEAD)', dset);
        assert(exist(xmat,'file')>0, 'getTsBase:noXmat','missing design matrix %s', xmat);

        % baseline (polynomial) coefficient sub-bricks, identified by label
        lab        = strsplit(strtrim(afniSys(afni,['3dinfo -label ' dset])), '|');
        isBaseCoef = contains(lab,'Pol') & endsWith(lab,'_Coef');
        assert(any(isBaseCoef),'getTsBase:noBaseCoef','no Pol*_Coef baseline sub-bricks in %s', dset);
        baseBrick  = find(isBaseCoef) - 1;             % 0-based sub-brick indices
        baseName   = erase(lab(isBaseCoef),'_Coef');   % names to match against xmat columns

        % read those coef sub-bricks via a temp NIfTI -> [X Y Z nBase]
        tmp = [tempname '.nii.gz'];
        sel = ['[' regexprep(strtrim(num2str(baseBrick)),'\s+',',') ']'];
        afniSys(afni,['3dTcat -overwrite -prefix ' tmp ' ''' dset sel '''']);
        beta = niftiread(tmp); delete(tmp);

        % design matrix: numbers (via 1dcat) + header metadata (labels, run starts)
        Xfull    = str2num(afniSys(afni,['1dcat ''' xmat ''''])); %#ok<ST2NM>   % [nT x nReg]
        hdr      = fileread(xmat);
        tok      = regexp(hdr,'ColumnLabels\s*=\s*"([^"]*)"','tokens','once');
        colLab   = strtrim(strsplit(tok{1},';'));
        tok      = regexp(hdr,'RunStart\s*=\s*"([^"]*)"','tokens','once');
        runStart = str2double(strsplit(tok{1},','));   % 0-based row index of each run start
        nT       = size(Xfull,1);

        % baseline design columns, matched to the coef bricks by name (order preserved)
        [tf,loc] = ismember(baseName, colLab);
        assert(all(tf),'getTsBase:colMismatch','baseline coef labels missing from xmat ColumnLabels');
        Xbase = Xfull(:,loc);                          % [nT x nBase]

        % sanity: reconstructed length must match the ts (no censoring assumed)
        nTts = sum(cellfun(@(im) size(im,4), vessel.im.ts.im));
        assert(nT==nTts, 'getTsBase:lenMismatch', ...
            'design has %d rows but ts has %d frames (censored xmat?)', nT, nTts);
        assert(numel(runStart)==numel(vessel.im.ts.im), 'getTsBase:runMismatch', ...
            'design has %d runs but ts has %d', numel(runStart), numel(vessel.im.ts.im));

        % crop coefficients to the ROI, exactly as vessel.im.ts was cropped
        % (loader convention: permute(niftiread,[2 1 3 4]) then index (y, x))
        xr = vessel.im.ts.x; yr = vessel.im.ts.y;
        beta = permute(beta,[2 1 3 4]);
        beta = beta(yr(1):yr(2), xr(1):xr(2), 1, :);
        ny = yr(2)-yr(1)+1; nx = xr(2)-xr(1)+1;
        beta = reshape(beta, ny*nx, []);               % [nVox x nBase]

        % full baseline timeseries (per voxel), split per run
        baseFull = beta * Xbase.';                      % [nVox x nT]
        runEnd   = [runStart(2:end) nT];                % 1-based slice ends
        tsBase       = vessel.im.ts;                    % mirror the ts struct layout
        tsBase.fName = {[dset '.BRIK']; [dset '.HEAD']; xmat};
        tsBase.im    = cell(1,numel(runStart));
        for r = 1:numel(runStart)
            cols = runStart(r)+1 : runEnd(r);
            tsBase.im{r} = reshape(baseFull(:,cols), ny, nx, 1, numel(cols));
        end
    end

    function out = afniSys(afni, cmd)
        [s,out] = system([afni '; ' cmd]);
        if s; error('getTsBase:afni','AFNI command failed:\n  %s\n%s', cmd, out); end
    end
end
