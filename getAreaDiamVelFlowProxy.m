function [vessel,fAll] = getAreaDiamVelFlowProxy(vessel,lumenMask,surroundMask,srcField,tValAvFlag)
    % Extract area, diameter, velocity and volumetric-flow proxy timeseries from
    % vessel images, for a caller-specified source field (the raw timeseries
    % vessel.im.ts, or a deconvolved response vessel.im.resp / vessel.im.resp2),
    % reusing a single core computation.
    %
    % srcField : source field name, or cellstr of field names, to process
    %            (default {'ts'}). For each field <srcField> this
    %            produces the absolute proxies
    %   vessel(v).im.<srcField>Area : area     proxy
    %   vessel(v).im.<srcField>Vel  : velocity proxy
    %   vessel(v).im.<srcField>D    : diameter proxy, D = 2*sqrt(A/pi)
    %                                 (circular cross-section assumption)
    %   and the fractional-change / volumetric-flow proxies (dX/X), all relative
    %   to a per-signal baseline X0:
    %   vessel(v).im.<srcField>AoA  : dA/A = (A-A0)/A0
    %   vessel(v).im.<srcField>VoV  : dV/V = (V-V0)/V0
    %   vessel(v).im.<srcField>DoD  : dD/D = (D-D0)/D0
    %   vessel(v).im.<srcField>QoQe : dQ/Q exact      = (1+dV/V).*(1+dA/A) - 1   (Q = V*A)
    %   vessel(v).im.<srcField>QoQa : dQ/Q 1st-order  = dV/V + dA/A
    %   Baseline X0: the proxy evaluated on the dedicated baseline image
    %   (vessel.im.basePolyRun, the .vecBase field) -- run-averaged for the
    %   deconvolved responses (resp*), per-run for the raw timeseries (ts). See
    %   dQoQ_derivation.md for the dQ/Q derivation.
    %
    % lumenMask / surroundMask : vessel.polyLabel names selecting the intravascular
    %            ("lumen", default 'peakVox') and surrounding ("surround", default
    %            'dilate1p5') voxel masks used to form the proxies. Valid polyLabel
    %            entries (built by modifyRoi; depends on the steps applied):
    %            'original', 'peakVox', 'dilate1', 'dilate1p5', 'dilate2', 'tissue'
    %            ('tissue' is also used internally as the tissue baseline).
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
    if nargin<2 || isempty(lumenMask);    lumenMask    = 'peakVox';   end
    if nargin<3 || isempty(surroundMask); surroundMask = 'dilate1p5'; end
    if nargin<4 || isempty(srcField);     srcField     = {'ts'};      end
    if nargin<5; tValAvFlag = true; end
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
             vessel(v).im.([fld 'Diam']),vessel(v).im.([fld 'Area']).wMask,vessel(v).im.([fld 'Area']).zMask,vessel(v).im.([fld 'Area']).tMask] = doIt(vessel(v),lumenMask,surroundMask,fld,addBaseFlag);

            % fold the fractional changes (dX/X) and volumetric flow change (dQ/Q)
            % into the proxy output. Baseline X0 = the proxy on the dedicated
            % baseline image (.vecBase). Alternative baselines (first frame /
            % per-run temporal mean) are available via baseMode in getdXoX.
            baseMode = 'vecBase';% 'mean'; if addBaseFlag; baseMode = 'first'; end
            vessel(v).im.([fld 'AoA' ])  = getdXoX(vessel(v).im.([fld 'Area']), baseMode, 'dA/A = (A-A0)/A0');
            vessel(v).im.([fld 'VoV' ])  = getdXoX(vessel(v).im.([fld 'Vel' ]), baseMode, 'dV/V = (V-V0)/V0');
            vessel(v).im.([fld 'DoD' ])  = getdXoX(vessel(v).im.([fld 'Diam']), baseMode, 'dD/D = (D-D0)/D0');

            vessel(v).im.([fld 'QoQ']) = getdQoQ2(vessel(v).im.([fld 'VoV' ]),vessel(v).im.([fld 'DoD']));
        end
    end

    function [outArea,outVel,outDiam,wMask,zMask,tMask] = doIt(vessel,lumenMask,surroundMask,fld,addBaseFlag)
        if ~exist('lumenMask'   ,'var') || isempty(lumenMask   ); lumenMask    = 'peakVox'; end
        if ~exist('surroundMask','var') || isempty(surroundMask); surroundMask = 'dilate1p5'; end
        wMask        = vessel.polyMask{ismember(vessel.polyLabel,lumenMask   )};
        zMask        = vessel.polyMask{ismember(vessel.polyLabel,surroundMask)};
        zMask(wMask) = false;
        tMask        = vessel.polyMask{ismember(vessel.polyLabel,'tissue'   )};

        %-----
        % since the surround mask is dilated from the lumen mask, it may contain tissue voxels
        % (the latter is computed separately). Let's just remove those voxels from the tissue mask.
        tMask(zMask) = false;
        %-----

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
        % outArea.wMask = wMask;
        % outArea.zMask = zMask;
        % outArea.tMask = tMask;
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
        % outVel.wMask = wMask;
        % outVel.zMask = zMask;
        % outVel.tMask = tMask;
        outVel.im = [];
        outVel.im2vec = wMask;
        outVel.info = 'vox x time';
        outVel.info2 = 'Sw: mean over intravascular voxels (actually just the peak voxel for now)';

        outDiam = vessel.im.(fld);
        outDiam.fName = '';
        % outDiam.wMask = wMask;
        % outDiam.zMask = zMask;
        % outDiam.tMask = tMask;
        outDiam.im = [];
        outDiam.im2vec = [];
        outDiam.info = 'vox x time';
        outDiam.info2 = 'D = 2*sqrt(A/pi): vessel diameter from area, assuming a circular cross-section';

        area = cell(1,nRun); areaBase = cell(1,nRun);
        vel  = cell(1,nRun); velBase  = cell(1,nRun);
        diam = cell(1,nRun); diamBase = cell(1,nRun);
        % frac = cell(1,nRun);
        base = vessel.im.basePolyRun.im;
        if strcmp(fld,'resp')
            % By the FIR deconvolution design, responses are assumed on the same scale across runs.
            % Appropriately using individual-run baseline would require scaling based on a proper reference (surrounding tissue?).
            [base{:}] = deal(mean(cat(3,base{:}),3));
        elseif strcmp(fld,'ts')
        else; error('getAreaDiamVelFlowProxy:badField','fld must be ''resp'' or ''ts''');
        end
        
        ~strcmp(fld,'resp') && ~strcmp(fld,'ts') && error('getAreaDiamVelFlowProxy:badField','fld must be ''resp'' or ''ts''');
        for r = 1:nRun
            if addBaseFlag
                im{r} = im{r} + base{r};
            end
            [area{r}    , vel{r}    , diam{r}    , frac] = computeAVD(im{r}  ,wMask,zMask,tMask,wN,zN,tValAvFlag,areaMethod);
            [areaBase{r}, velBase{r}, diamBase{r}, frac] = computeAVD(base{r},wMask,zMask,tMask,wN,zN,tValAvFlag,areaMethod);

            if any(frac>1) || any(frac<0)
                warning('getAreaDiamVelFlowProxy:fOutOfBounds', ...
                    '%s, vessel %d, run %d: %d/%d surround-fraction samples outside [0,1] (%d<0, %d>1).', ...
                    fld,v,r,nnz(frac<0 | frac>1),numel(frac),nnz(frac<0),nnz(frac>1));
            end
        end

        % single-image sources collapse the per-run cell back to a matrix
        if ~isCell
            outArea.vec     = area{1};     outVel.vec     = vel{1};     outDiam.vec     = diam{1};
            outArea.vecBase = areaBase{1}; outVel.vecBase = velBase{1}; outDiam.vecBase = diamBase{1};
        else
            outArea.vec     = area;     outVel.vec     = vel;     outDiam.vec     = diam;
            outArea.vecBase = areaBase; outVel.vecBase = velBase; outDiam.vecBase = diamBase;
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
        frac = min(max(frac,0),1);
        % frac(frac<0 | frac>1) = nan;
        % frac(frac<0) = nan;
        switch areaMethod
            case 1 % Variant 1
                area = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);
            case 2 % Variant 2: allows boundind the surround fraction (zVal-tVal)./(wVal-tVal) to [0,1]
                area = wN + zN.*frac;
            otherwise
                error('getAreaDiamVelFlowProxy:badAreaMethod','areaMethod must be 1 or 2');
        end
        area = permute(area,[2 1 3 4]);
        vel  = permute(wVal,[2 1 3 4]);
        diam = 2*sqrt(area/pi);
    end

    % --- fractional change dX/X of a proxy's .vec, preserving struct metadata.
    %     baseMode 'first' -> baseline = first time frame; 'mean' -> per-run temporal
    %     mean; 'vecBase' -> the proxy's dedicated baseline (.vecBase, per element).
    function s = getdXoX(proxy,baseMode,info2)
        s = proxy; s.info2 = info2;
        switch baseMode
            case 'first';   base = @(x) x(:,1);
            case 'mean';    base = @(x) mean(x,2,'omitnan');
            case 'vecBase'; base = @(x) x;
            otherwise; error('getdXoX:badBaseMode','baseMode must be ''first'', ''mean'' or ''vecBase''');
        end
        dXoX = @(x,b) (x - base(b))./base(b);
        if iscell(proxy.vec)
            s.vec = cellfun(dXoX,proxy.vec,proxy.vecBase,'UniformOutput',false);
        else
            s.vec = dXoX(proxy.vec,proxy.vecBase);
        end
    end

    function QoQ = getdQoQ2(VoV,DoD)
        QoQ = VoV;
        for r = 1:length(VoV.vec)
            QoQ.vec{r}     = computeQoQ(VoV.vec{r}    ,DoD.vec{r}    );
            QoQ.vecBase{r} = [];
        end
        QoQ.info2 = 'QoQ = (1+VoV).*(1+DoD).^2 - 1;';
    end


end
