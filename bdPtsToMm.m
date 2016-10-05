function bdPtsMm = bdPtsToMm(bd_pts,patchSizeMm)

% function bdPtsToMm 
% This function converts the boundary points from Aura Tools (in units of
% pixels) to mm by scaling by the width of the B-Scan.
%
% Inputs:
% bd_pts      = The matrix of boundary points outputed from Aura Tools
% patchSizeMm = The width of the B-Scan in mm. DEFAULT = 6mm which is what  
%               Aura Tools cuts the OCT data to. 
%
% MAB OCT 2016

if isempty(patchSizeMm)
    patchSizeMm = 6;
end

numPxls = size(bd_pts,1);
pxlSize = patchSizeMm/numPxls;
bdPtsMm = bd_pts.*pxlSize;

end
