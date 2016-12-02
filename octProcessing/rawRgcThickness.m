function rgcPlus = rawRgcThickness(bd_pts,header)
% This function takes in the segmentation from Aura Tools and uses the 
% header to compute the depth of each boundary layer in mm.
% The output is the subtraction of NFL boundary the RGC/IPL boundary 
% resulting in the RGC+ layer thickness. 


bdPtsMm = bdPtsToMm(bd_pts,header); % convert pixels to mm
rgcPlus = bdPtsMm(:,:,3)-bdPtsMm(:,:,2); % Get the RGC+ layer thickness

end

