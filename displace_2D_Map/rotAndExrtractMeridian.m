function [RGCdenisty_mmSq RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDeg)

%rotate and crop
RF_crop   = imrotate(RFdensity,rotDeg,'crop','bilinear');
rotRGC_crop  = imrotate(RGCdensity,rotDeg,'crop','bilinear');

%% Get meridian fro RFdensity
RF_midPoint = round(

end