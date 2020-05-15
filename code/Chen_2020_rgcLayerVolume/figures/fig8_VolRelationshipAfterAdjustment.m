function fig8_VolRelationshipAfterAdjustment(XPos_Degs,comboTable,scoreExpandedSmoothed,adjustedCoeff,coeff,gcVolumePerDegSq,nDimsToUse,saveDir)

% Plot the synthesized reconstructions by axial length
gcVec_reconstruct = scoreExpandedSmoothed(:,1:nDimsToUse)*adjustedCoeff(:,1:nDimsToUse)';
gcVec_reconstruct(isnan(gcVolumePerDegSq))=nan;
gcVec_reconstruct_orig = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff(:,1:nDimsToUse)';
meangcVec_reconstruct = nanmean(gcVec_reconstruct,2);
gcMedianVol_reconstruct = nanmedian(gcVec_reconstruct,1);
gcMeanVol_reconstruct = nanmean(gcVec_reconstruct,1);
gcMeanVol_reconstruct_orig = nanmean(gcVec_reconstruct_orig,1);

% Censor any reconstructed value that is less than zero
gcVec_reconstruct(gcVec_reconstruct<0)=nan;

% Plot the GC thickness functions
str = ['Adjusted GC vol profiles for each subject (and mean)'];
h=profilePlot(XPos_Degs, gcVec_reconstruct, meangcVec_reconstruct, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig8','a.png'));

% Plot GC thickness vs axial length, 
str = ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,gcMedianVol_reconstruct'))];
h=regressionPlot(comboTable.Axial_Length_average, gcMedianVol_reconstruct', 'Axial Length [mm]','Median GC Tissue Volume [mm^3 / deg^2]', [],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig8','b.png'));

% Plot GC thickness vs axial length
str = ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,gcMeanVol_reconstruct'))];
h=regressionPlot(comboTable.Axial_Length_average, gcMeanVol_reconstruct', 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', [],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig8','c.png'));

% Plot GC thickness vs axial length
str = ['Axial length vs. original mean GC vol, r=',num2str(corr(comboTable.Axial_Length_average,gcMeanVol_reconstruct_orig'))];
h=regressionPlot(comboTable.Axial_Length_average, gcMeanVol_reconstruct_orig', 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', [],1);
title(str);
setTightFig
saveas(h,fullfile(saveDir,'fig8','d.png'));