function [vessel,fAll] = getAreaVel(vessel,lumenMask,surroundMask,srcField,tValAvFlag)
    % Extract area and velocity proxy timeseries from vessel images, for a
    % caller-specified source field (raw timeseries vessel.im.ts, or a deconvolved
    % response vessel.im.resp / resp2). Area+velocity-only sibling of
    % getAreaDiamVelFlowProxy: it shares the core area/velocity computation
    % (computeAreaVel, area variant 2 by default) but does NOT compute diameter, dD/D or dQ/Q
    % -- those are handled separately, later in the pipeline.
    %
    % srcField : source field name, or cellstr of field names (default {'ts'}). For
    %            each field <srcField> this produces the absolute proxies
    %   vessel(v).im.<srcField>Area : area     proxy (.vec and .vecBase)
    %   vessel(v).im.<srcField>Vel  : velocity proxy (.vec and .vecBase)
    %   The dedicated baseline image proxy is kept in .vecBase (per run for the raw
    %   timeseries ts, run-averaged for deconvolved responses resp) so fractional
    %   changes (dA/A, dV/V) can be formed later, where/when D, Q and Faa are.
    %
    % lumenMask / surroundMask : vessel.polyLabel names selecting the intravascular
    %            ("lumen", default 'peakVox') and surround ("surround", default
    %            'dilate1p5') voxel masks. See getAreaDiamVelFlowProxy for the valid
    %            polyLabel entries and the baseline-recovery / .vec layout details
    %            (mirrored here).
    if nargin<2 || isempty(lumenMask);    lumenMask    = 'peakVox';   end
    if nargin<3 || isempty(surroundMask); surroundMask = 'dilate1p5'; end
    if nargin<4 || isempty(srcField);     srcField     = {'ts'};      end
    if nargin<5; tValAvFlag = true; end
    if ~iscell(srcField); srcField = cellstr(srcField); end

    fAll = struct();
    for v = 1:length(vessel)
        for k = 1:numel(srcField)
            fld = srcField{k};
            if ~isfield(vessel(v).im,fld) || isempty(vessel(v).im.(fld)); continue; end
            % deconvolved responses need the baseline added back; raw ts does not
            addBaseFlag = startsWith(fld,'resp');
            [vessel(v).im.([fld 'Area']), ...
             vessel(v).im.([fld 'Vel']), ...
             vessel(v).im.([fld 'Area']).wMask, ...
             vessel(v).im.([fld 'Area']).zMask, ...
             vessel(v).im.([fld 'Area']).tMask] = doIt(vessel(v),lumenMask,surroundMask,fld,addBaseFlag);
        end
    end

    function [outArea,outVel,wMask,zMask,tMask] = doIt(vessel,lumenMask,surroundMask,fld,addBaseFlag)
        wMask        = vessel.polyMask{ismember(vessel.polyLabel,lumenMask   )};
        zMask        = vessel.polyMask{ismember(vessel.polyLabel,surroundMask)};
        zMask(wMask) = false;
        tMask        = vessel.polyMask{ismember(vessel.polyLabel,'tissue'   )};
        % the surround mask is dilated from the lumen, so it may contain tissue
        % voxels (computed separately) -- remove those from the tissue mask.
        tMask(zMask) = false;

        % refactor the source images into a per-run cell (isCell flags the layout)
        im = vessel.im.(fld).im;
        isCell = iscell(im);
        if ~isCell; im = {im}; end
        nRun = numel(im);

        % shared output templates (carry over dt, masks, info, ...)
        outArea = vessel.im.(fld);
        outArea.fName = ''; outArea.im = []; outArea.im2vec = [];
        outArea.info  = 'vox x time';
        outArea.info2 = ['(Nw*(Sw-St)+Nz*(Sz-St)) / (Sw-St)' newline...
                         'Nw: number of intravascular voxels' newline...
                         'Nz: number of surrounding voxels' newline...
                         'Sw: mean signal in intravascular voxels' newline...
                         'Sz: mean signal in surrounding voxels' newline...
                         'St: mean signal in tissue voxels'];

        outVel = vessel.im.(fld);
        outVel.fName = ''; outVel.im = []; outVel.im2vec = wMask;
        outVel.info  = 'vox x time';
        outVel.info2 = 'Sw: mean over intravascular voxels (actually just the peak voxel for now)';

        area = cell(1,nRun); areaBase = cell(1,nRun);
        vel  = cell(1,nRun); velBase  = cell(1,nRun);
        base = vessel.im.basePolyRun.im;
        if strcmp(fld,'resp')
            % By the FIR deconvolution design, responses are assumed on the same
            % scale across runs -> use the run-average baseline.
            [base{:}] = deal(mean(cat(3,base{:}),3));
        elseif ~strcmp(fld,'ts')
            error('getAreaVel:badField','fld must be ''resp'' or ''ts''');
        end
        for r = 1:nRun
            if addBaseFlag; im{r} = im{r} + base{r}; end
            [area{r}    , vel{r}    , frac] = computeAreaVel(im{r}  ,wMask,zMask,tMask,tValAvFlag);
            [areaBase{r}, velBase{r}, ~   ] = computeAreaVel(base{r},wMask,zMask,tMask,tValAvFlag);

            if any(frac>1) || any(frac<0)
                warning('getAreaVel:fOutOfBounds', ...
                    '%s, vessel %d, run %d: %d/%d surround-fraction samples outside [0,1] (%d<0, %d>1).', ...
                    fld,v,r,nnz(frac<0 | frac>1),numel(frac),nnz(frac<0),nnz(frac>1));
            end
        end

        % single-image sources collapse the per-run cell back to a matrix
        if ~isCell
            outArea.vec = area{1}; outVel.vec = vel{1};
            outArea.vecBase = areaBase{1}; outVel.vecBase = velBase{1};
        else
            outArea.vec = area; outVel.vec = vel;
            outArea.vecBase = areaBase; outVel.vecBase = velBase;
        end
    end
end
