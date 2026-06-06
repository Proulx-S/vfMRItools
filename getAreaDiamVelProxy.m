function [vessel,fAll] = getAreaDiamVelProxy(vessel,srcField,tValAvFlag)
    % Extract area, diameter and velocity proxy timeseries from vessel images,
    % for a caller-specified source field (the raw timeseries vessel.im.ts, or a
    % deconvolved response vessel.im.resp / vessel.im.resp2), reusing a single
    % core computation.
    %
    % srcField : source field name, or cellstr of field names, to process
    %            (default {'ts'}). For each field <srcField> this
    %            produces
    %   vessel(v).im.<srcField>Area : area     proxy
    %   vessel(v).im.<srcField>Vel  : velocity proxy
    %   vessel(v).im.<srcField>D    : diameter proxy, D = 2*sqrt(A/pi)
    %                                 (circular cross-section assumption)
    %
    % Each source field may hold either a single image (e.g. an averaged
    % response: <srcField>.im is a numeric array) or one image per run (e.g. raw
    % runs or per-run responses: <srcField>.im is a cell). The output .vec
    % mirrors the source: a single [vox x time] matrix for a single image, or a
    % cell with one [vox x time] matrix per run otherwise.
    %
    % Deconvolved responses (resp, resp2) carry no absolute signal level, so the
    % per-run baseline (vessel.im.basePolyRun) is added back before forming the
    % proxies (mirroring getVesselResp). The raw timeseries (ts) already carries
    % absolute signal levels, so no baseline is added.
    %
    % NB: how the baseline is recovered (vessel.im.basePolyRun, with a per-run
    % image for per-run sources and the run-average for a single/averaged source)
    % is idiosyncratic to our specific data structure -- see getBase below.
    if nargin<2 || isempty(srcField); srcField = {'ts'}; end
    if nargin<3; tValAvFlag = true; end
    if ~iscell(srcField); srcField = cellstr(srcField); end

    % Area proxy method:
    %   1: ( Nw*(Sw-St) + Nz*(Sz-St) ) / (Sw-St)
    %   2: Nw + Nz*frac, with frac = (Sz-St)/(Sw-St) bounded to [0,1]
    areaMethod = 2;

    fAll = struct();
    for v = 1:length(vessel)
        for k = 1:numel(srcField)
            fld = srcField{k};
            if ~isfield(vessel(v).im,fld) || isempty(vessel(v).im.(fld)); continue; end
            % deconvolved responses need the baseline added back; raw ts does not
            addBaseFlag = startsWith(fld,'resp');
            [vessel(v).im.([fld 'Area']), ...
             vessel(v).im.([fld 'Vel']), ...
             vessel(v).im.([fld 'D']), ...
             fAll(v).(fld)] = doIt(vessel(v),fld,addBaseFlag);
        end
    end

    function [outArea,outVel,outDiam,frac] = doIt(vessel,fld,addBaseFlag)
        wMask        = vessel.polyMask{ismember(vessel.polyLabel,'peakVox'  )};
        zMask        = vessel.polyMask{ismember(vessel.polyLabel,'dilate1p5')};
        zMask(wMask) = false;
        tMask        = vessel.polyMask{ismember(vessel.polyLabel,'tissue'   )};
        wN = nnz(wMask);
        zN = nnz(zMask);

        % refactor the source images into a per-run cell (isCell flags the layout)
        im = vessel.im.(fld).im;
        isCell = iscell(im);
        if ~isCell; im = {im}; end
        nRun = numel(im);

        % shared output templates (carry over dt, masks, info, ...)
        outArea = vessel.im.(fld);
        outArea.fName = '';
        outArea.maskResp.wMask = wMask;
        outArea.maskResp.zMask = zMask;
        outArea.maskResp.tMask = tMask;
        outArea.im = [];
        outArea.im2vec = [];
        outArea.info = 'vox x time';
        outArea.info2 = ['(Nw*(Sw-St)+Nz*(Sz-St)) / (Sw-St)' newline...
                                   'Nw: number of intravascular voxels' newline...
                                   'Nz: number of surrounding voxels' newline...
                                   'Sw: mean signal in intravascular voxels' newline...
                                   'Sz: mean signal in surrounding voxels' newline...
                                   'St: mean signal in tissue voxels'   ];

        outVel = vessel.im.(fld);
        outVel.fName = '';
        outVel.maskResp.wMask = wMask;
        outVel.maskResp.zMask = zMask;
        outVel.maskResp.tMask = tMask;
        outVel.im = [];
        outVel.im2vec = wMask;
        outVel.info = 'vox x time';
        outVel.info2 = 'Sw: mean over intravascular voxels (actually just the peak voxel for now)';

        outDiam = vessel.im.(fld);
        outDiam.fName = '';
        outDiam.maskResp.wMask = wMask;
        outDiam.maskResp.zMask = zMask;
        outDiam.maskResp.tMask = tMask;
        outDiam.im = [];
        outDiam.im2vec = [];
        outDiam.info = 'vox x time';
        outDiam.info2 = 'D = 2*sqrt(A/pi): vessel diameter from area, assuming a circular cross-section';

        area = cell(1,nRun);
        vel  = cell(1,nRun);
        diam = cell(1,nRun);
        frac = cell(1,nRun);
        for r = 1:nRun
            base = getBase(vessel,addBaseFlag,isCell,r);
            [area{r},vel{r},diam{r},frac{r}] = computeAVD(im{r}+base,wMask,zMask,tMask,wN,zN,tValAvFlag,areaMethod);
            if any(frac{r}>1) || any(frac{r}<0)
                warning('getAreaDiamVelProxy:fOutOfBounds', ...
                    '%s, vessel %d, run %d: %d/%d surround-fraction samples outside [0,1] (%d<0, %d>1).', ...
                    fld,v,r,nnz(frac{r}<0 | frac{r}>1),numel(frac{r}),nnz(frac{r}<0),nnz(frac{r}>1));
            end
        end

        % single-image sources collapse the per-run cell back to a matrix
        if ~isCell
            outArea.vec = area{1}; outVel.vec = vel{1}; outDiam.vec = diam{1};
            frac = frac{1};
        else
            outArea.vec = area; outVel.vec = vel; outDiam.vec = diam;
        end
    end

    % per-run baseline to add back to a deconvolved response (0 for raw ts)
    function base = getBase(vessel,addBaseFlag,isCell,r)
        if ~addBaseFlag; base = 0; return; end
        if isCell
            % per-run source -> per-run baseline
            base = vessel.im.basePolyRun.im{r};
        else
            % single/averaged source -> baseline averaged across runs
            base = mean(cat(3,vessel.im.basePolyRun.im{:}),3);
        end
    end

    % --- shared core: the Area, Vel and D computation, reused by every scenario
    function [area,vel,diam,frac] = computeAVD(im4d,wMask,zMask,tMask,wN,zN,tValAvFlag,areaMethod)
        tsIm = permute(im4d,[4 1 2 3]); % time x X x Y x Z
        wVal = mean(tsIm(:,wMask),2);
        zVal = mean(tsIm(:,zMask),2);
        tVal = mean(tsIm(:,tMask),2);
        if tValAvFlag
            tVal = mean(tVal,1); % assume stable tissue signal to avoid noise
        end
        % surround fraction (unbounded; reported for diagnostics)
        frac = (zVal-tVal)./(wVal-tVal);
        switch areaMethod
            case 1 % Variant 1
                area = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);
            case 2 % Variant 2: bound the surround fraction to [0,1]
                area = wN + zN.*min(max(frac,0),1);
            otherwise
                error('getAreaDiamVelProxy:badAreaMethod','areaMethod must be 1 or 2');
        end
        area = permute(area,[2 1 3 4]);
        vel  = permute(wVal,[2 1 3 4]);
        diam = 2*sqrt(area/pi);
    end
end
