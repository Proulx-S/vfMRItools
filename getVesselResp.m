function vessel = getVesselResp(vessel)
    % Uresp: temporal singular vectors scaled by singular values [copmonents x time]
    % Vresp: spatial singular vectors scaled by singular values  [components x X x Y]
    % AreaResp: area response     [1 x time]
    % VelResp : velocity response [1 x time]
    % PeakVoxResp: peak voxel response     [1   x time]
    % SurrVoxResp: surround voxel response [vox x time]

    for v = 1:length(vessel)
        [vessel(v).svdResp,vessel(v).im.respArea,vessel(v).im.respVel,vessel(v).im.respPeakVox,vessel(v).im.respSurVox] = doIt(vessel(v));
    end

    function [svdResp,respArea,respVel,respPeakVox,respSurVox] = doIt(vessel)
        respIm = permute(vessel.im.resp.im,[4 1 2 3]);
        wMask = vessel.polyMask{ismember(vessel.polyLabel,'peakVox')};
        zMask = vessel.polyMask{ismember(vessel.polyLabel,'dilate1p5')}; zMask(wMask) = false;
        tMask = vessel.polyMask{ismember(vessel.polyLabel,'tissue')};
        d2Mask = vessel.polyMask{ismember(vessel.polyLabel,'dilate2')};
        
        % SVD transform
        svdMask = tMask|d2Mask;
        [U,S,V] = svd(respIm(:,svdMask),'econ','vector'); % time x vox (excluding those containing other vessels)
        % rectify singular vectors for a positive maximum deflection of the temporal singular vector
        [~,b] = max(abs(U),[],1);
        flp = sign(U(sub2ind(size(U),b,1:size(U,2))));
        U = U.*flp;
        V = V.*flp;
        % output
        svdResp = vessel.im.base;
        svdResp.fName   = '';
        svdResp.im      = [];
        svdResp.maskSVD = svdMask;
        svdResp.sv      = S;
        svdResp.svSpace = permute(V,[2 1]).*S;
        svdResp.svTime  = permute(U,[2 1]).*S;
        svdResp.info   = 'component x vox/time';
        svdResp.info2 = ['sv:      singular values' newline 'svSpace: spatial singular vectors scaled by singular values' newline 'svTime:  temporal singular vectors scaled by singular values'];
        
        % Area/velocity transform
        base = permute(mean(cat(3,vessel.im.basePolyRun.im{:}),3),[4 1 2 3]);
        respIm = respIm + base;
        wVal = mean(respIm(:,wMask),2); wN = nnz(wMask);
        zVal = mean(respIm(:,zMask),2); zN = nnz(zMask);
        tVal = mean(respIm(:,tMask),2); tN = nnz(tMask);
        tVal = mean(tVal,1); % assume stable tissue signal to avoid noise
        f = (zVal-tVal)./(wVal-tVal);
        f = min(max(f,0),1); % bound f from 0 to 1.
        AreaResp = wN + zN.*f;    
        % AreaResp = ( wN.*(wVal-tVal) + zN.*(zVal-tVal) ) ./ (wVal-tVal);
        respArea = vessel.im.resp;
        respArea.fName = '';
        respArea.maskResp.wMask = wMask;
        respArea.maskResp.zMask = zMask;
        respArea.maskResp.tMask = tMask;
        respArea.im = [];
        respArea.im2vec = [];
        respArea.vec = permute(AreaResp,[2 1]);
        respArea.info = 'vox x time';
        respArea.info2 = ['(Nw*(Sw-St)+Nz*(Sz-St)) / (Sw-St)' newline...
                                   'Nw: number of intravascular voxels' newline...
                                   'Nz: number of surrounding voxels' newline...
                                   'Sw: mean signal in intravascular voxels' newline...
                                   'Sz: mean signal in surrounding voxels' newline...
                                   'St: mean signal in tissue voxels'   ];
        
        respVel = vessel.im.resp;
        respVel.fName = '';
        respVel.maskResp.wMask = wMask;
        respVel.maskResp.zMask = zMask;
        respVel.maskResp.tMask = tMask;
        respVel.im = [];
        respVel.im2vec = wMask;
        respVel.vec = permute(mean(wVal,2),[2 1 3 4]);
        respVel.info = 'vox x time';
        respVel.info2 = 'Sw: mean over intravascular voxels (actually just the peak voxel for now)';
        
        % Peak/surround voxel
        respPeakVox = vessel.im.resp;
        respPeakVox.fName = '';
        respPeakVox.maskResp = vessel.polyMask{ismember(vessel.polyLabel,'peakVox')};
        respPeakVox.im = [];
        respPeakVox.im2vec = respPeakVox.maskResp;
        respPeakVox.vec = permute(respIm(:,respPeakVox.im2vec),[2 1]);
        respPeakVox.info = 'vox x time';
        respPeakVox.info2 = 'peak signal intravascular voxel';
        
        respSurVox = vessel.im.resp;
        respSurVox.fName = '';
        sMask = vessel.polyMask{ismember(vessel.polyLabel,'dilate1')};
        sMask(vessel.polyMask{ismember(vessel.polyLabel,'original')}) = false;
        respSurVox.maskResp = sMask;
        respSurVox.im = [];
        respSurVox.im2vec = respSurVox.maskResp;
        respSurVox.vec = permute(respIm(:,respSurVox.im2vec),[2 1]);
        respSurVox.info = 'vox x time';
        respSurVox.info2 = 'surround voxels';