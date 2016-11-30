function [RGCdenisty_mmSq,RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDeg)

%% rotate and crop
RF_crop   = imrotate(RFdensity,-1.*rotDeg,'crop','bilinear');
RGC_crop  = imrotate(RGCdensity,-1.*rotDeg,'crop','bilinear');

%% Get meridian from RF density rotated input image
RF_midPoint = round(size(RF_crop,1)/2);
RFdensity_sqDeg = RF_crop(RF_midPoint,RF_midPoint:end);
%% Get meridian from RGC density rotated input image
RGC_midPoint = round(size(RGC_crop,1)/2);
RGCdenisty_mmSq = RF_crop(RGC_midPoint,RGC_midPoint:end);

end