function [area,vel,frac] = computeAreaVel(im4d,wMask,zMask,tMask,tValAvFlag)
    % Area and velocity proxies for one vessel image, extracted from the shared
    % core of getAreaDiamVelFlowProxy (computeAVD), restricted to area & velocity
    % (no diameter) and using area variant 1.
    %
    %   im4d       : X x Y x Z x time image
    %   wMask      : intravascular ("lumen") voxel mask
    %   zMask      : surrounding voxel mask (lumen voxels already removed)
    %   tMask      : tissue voxel mask (surround voxels already removed)
    %   tValAvFlag : true (default) -> average the tissue signal over time
    %                (stable-tissue assumption, reduces noise); false -> keep it
    %                time-resolved
    %
    %   area : 1 x time area proxy, variant 1:
    %          area = ( Nw*(Sw-St) + Nz*(Sz-St) ) / (Sw-St)
    %   vel  : 1 x time velocity proxy = Sw (mean lumen signal)
    %   frac : 1 x time surround fraction (Sz-St)/(Sw-St), unbounded (diagnostic;
    %          variant 1 does not use it)
    %
    %   Nw, Nz are the voxel counts of wMask, zMask (computed here). Sw, Sz, St are
    %   the mean signals over the lumen, surround and tissue voxels.
    if ~exist('tValAvFlag','var') || isempty(tValAvFlag); tValAvFlag = true; end
    wN = nnz(wMask);
    zN = nnz(zMask);
    tsIm = permute(im4d,[4 1 2 3]);     % time x X x Y x Z
    wVal = mean(tsIm(:,wMask),2);
    zVal = mean(tsIm(:,zMask),2);
    tVal = mean(tsIm(:,tMask),2);
    if tValAvFlag
        tVal = mean(tVal,1);            % assume stable tissue signal to avoid noise
    end
    frac = (zVal-tVal)./(wVal-tVal);    % surround fraction (unbounded; diagnostic)
    area = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);   % variant 1
    area = permute(area,[2 1 3 4]);
    vel  = permute(wVal,[2 1 3 4]);
end
