function [rotDf_crop,rotRGC_crop] = rotAndCrop(Df,dRGC,rotDeg)


%rotate and crop
rotDf_crop   = imrotate(Df,rotDeg,'crop','bilinear');
rotRGC_crop  = imrotate(dRGC,rotDeg,'crop','bilinear');

end