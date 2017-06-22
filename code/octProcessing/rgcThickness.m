function [rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header)
% This function takes in the segmentation from Aura Tools and uses the 
% header to compute the depth of each boundary layer in mm.
% The output is the subtraction of NFL boundary the RGC/IPL boundary 
% resulting in the RGC+ layer thickness. 


bdPtsMm = bdPtsToMm(bd_pts,header); % convert pixels to mm
rgcPlusRaw = bdPtsMm(:,:,3)-bdPtsMm(:,:,2); % Get the RGC+ layer thickness

%% resample the RGC+ map to matche the resolution of the number of b-scans per scan

rgcPlus = imresize(rgcPlusRaw,[size(rgcPlusRaw,2),size(rgcPlusRaw,2)],'nearest');

%% Create Smaple Base
numSmps = size(rgcPlus,1);
image_size = header.SizeX * header.ScaleX;
smps = image_size./numSmps;
sampleBaseRadius = 0:smps:(image_size/2);

end

