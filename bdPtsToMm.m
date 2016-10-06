function bdPtsMm = bdPtsToMm(bd_pts,header)

% function bdPtsToMm 
% This function converts the boundary points from Aura Tools (in units of
% pixels) to mm by scaling by the width of the B-Scan.
%
% Inputs:
% bd_pts  = The matrix of boundary points outputed from Aura Tools
% header  = The accociated header file that is from the output of Arua
%           Tools. This is used to get the mm size of each pixel.
%
% MAB OCT 2016

if isempty(patchSizeMm)
    patchSizeMm = 6;
end

pxlSize = header.ScaleX;
bdPtsMm = bd_pts.*pxlSize;

end
