function rgcPlus = rawRgcThickness(bd_pts,header)
% This function takes in the output 

bdPtsMm = bdPtsToMm(bd_pts,header)
rgclayer = bdPtsMm(:,:,3)-bdPtsMm(:,:,2);

