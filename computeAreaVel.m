function [area,vel,frac] = computeAreaVel(im4d,wMask,zMask,tMask,tValAvFlag,areaMethod)
    % Area and velocity proxies for one vessel image, extracted from the shared
    % core of getAreaDiamVelFlowProxy (computeAVD), restricted to area & velocity
    % (no diameter).
    %
    %   im4d       : X x Y x Z x time image
    %   wMask      : intravascular ("lumen") voxel mask
    %   zMask      : surrounding voxel mask (lumen voxels already removed)
    %   tMask      : tissue voxel mask (surround voxels already removed)
    %   tValAvFlag : true (default) -> average the tissue signal over time
    %                (stable-tissue assumption, reduces noise); false -> keep it
    %                time-resolved
    %   areaMethod : area proxy variant (default 2)
    %                1: ( Nw*(Sw-St) + Nz*(Sz-St) ) / (Sw-St)
    %                2: Nw + Nz*frac, with frac bounded to [0,1]
    %
    %   area : 1 x time area proxy (per areaMethod)
    %   vel  : 1 x time velocity proxy = Sw (mean lumen signal)
    %   frac : 1 x time surround fraction (Sz-St)/(Sw-St), UNbounded (diagnostic)
    %
    %   Nw, Nz are the voxel counts of wMask, zMask (computed here). Sw, Sz, St are
    %   the mean signals over the lumen, surround and tissue voxels.
    if ~exist('tValAvFlag','var') || isempty(tValAvFlag); tValAvFlag = true; end
    if ~exist('areaMethod','var') || isempty(areaMethod); areaMethod = 2;    end
    wN = nnz(wMask);
    zN = nnz(zMask);
    tsIm = permute(im4d,[4 1 2 3]);     % time x X x Y x Z
    wVal = mean(tsIm(:,wMask),2);
    zVal = mean(tsIm(:,zMask),2);
    tVal = mean(tsIm(:,tMask),2);
    if tValAvFlag
        tVal = mean(tVal,1,'omitnan'); % assume stable tissue signal to avoid noise (omitnan: ignore censored frames)
    end
    frac = (zVal-tVal)./(wVal-tVal);    % surround fraction (unbounded; diagnostic)
    switch areaMethod
        case 1
            area = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);
        case 2
            area = wN + zN.*min(max(frac,0),1);          % bound the surround fraction to [0,1]
        otherwise
            error('computeAreaVel:badAreaMethod','areaMethod must be 1 or 2');
    end
    area = permute(area,[2 1 3 4]);
    vel  = permute(wVal,[2 1 3 4]);
end
