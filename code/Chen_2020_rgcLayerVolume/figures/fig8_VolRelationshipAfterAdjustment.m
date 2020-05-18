function fig8_VolRelationshipAfterAdjustment(XPos_Degs,comboTable,scoreExpandedSmoothed,adjustedCoeff,coeff,gcVolumePerDegSq,nDimsToUse,saveDir)

% Plot the synthesized reconstructions by axial length
gcVec_reconstruct = scoreExpandedSmoothed(:,1:nDimsToUse)*adjustedCoeff(:,1:nDimsToUse)';
meangcVolVec_reconstruct = nanmean(gcVec_reconstruct,2);
gcMeanVol_reconstruct = nanmean(gcVec_reconstruct,1);

% Censor any reconstructed value that is less than zero
gcVec_reconstruct(gcVec_reconstruct<0)=nan;

% Plot the GC volume functions, adjusted for axial length
str = ['Adjusted GC vol profiles for each subject (and mean)'];
h=profilePlot(XPos_Degs, gcVec_reconstruct, meangcVolVec_reconstruct, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig8','a.pdf'));

% Plot adjusted GC vol vs axial length
str = ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,gcMeanVol_reconstruct'))];
h=regressionPlot(comboTable.Axial_Length_average, gcMeanVol_reconstruct', 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', [],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig8','c.pdf'));
