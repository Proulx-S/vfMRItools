function cnsr = getCnsr(fList, imTemplate)
    % Per-run censoring info, packaged in the vessel.im.ts format but with the data
    % carried in .vec (per-run row vectors), like the other proxy .vec fields.
    %
    % fList : per-run cell of censor .csv paths (e.g. rCond.fCnsr or
    %         rCond.fCnsr_mainClust). Each csv is [nFrame x 2]: col 1 = frame index,
    %         col 2 = INCLUSION flag in the file (0 = exclude, 1 = include).
    % imTemplate : a vessel.im.ts struct whose layout / per-run frame counts is mirrored.
    %
    % Output cnsr : like vessel.im.ts but with the data in .vec instead of .im
    %   .fName : per-run censor csv paths (fList)
    %   .x/.y/.mask : copied from imTemplate (grid/run context)
    %   .vec : per-run cell, each [1 x nFrame] logical CENSOR mask -- true = EXCLUDE.
    %          The file's inclusion flag is inverted so the value matches the 'cnsr'
    %          name: 1 = censored/excluded, 0 = kept.

    cnsr       = imTemplate;
    cnsr.fName = fList;
    if isfield(cnsr,'im'); cnsr = rmfield(cnsr,'im'); end
    cnsr.vec   = cell(size(imTemplate.im));
    for r = 1:numel(fList)
        c  = readmatrix(fList{r});
        nT = size(imTemplate.im{r},4);
        assert(size(c,1)==nT, 'getCnsr:lenMismatch', ...
            'censor file %s has %d frames but ts run %d has %d.', fList{r}, size(c,1), r, nT);
        cnsr.vec{r} = reshape(~logical(c(:,2)),1,nT);   % invert: 1 = censored/excluded, 0 = kept
    end
end
