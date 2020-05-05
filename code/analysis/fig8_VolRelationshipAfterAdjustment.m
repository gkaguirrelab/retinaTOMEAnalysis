% Plot the synthesized reconstructions by axial length
gcVec_reconstruct = scoreExpandedSmoothed(:,1:nDimsToUse)*adjustedCoeff(:,1:nDimsToUse)';
gcVec_reconstruct_orig = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff(:,1:nDimsToUse)';
meangcVec_reconstruct = mean(gcVec_reconstruct,2);
gcMedianThick_reconstruct = nanmedian(gcVec_reconstruct,1);
gcMeanThick_reconstruct = nanmean(gcVec_reconstruct,1);
gcMeanThick_reconstruct_orig = nanmean(gcVec_reconstruct_orig,1);

% Plot the GC thickness functions, ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, gcVec_reconstruct, meangcVec_reconstruct, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig8','a.png'));


% Plot GC thickness vs axial length, ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick))]
h=regressionPlot(comboTable.Axial_Length_average, gcMedianThick_reconstruct', 'Axial Length [mm]','Median GC Tissue Volume [mm^3 / deg^2]', [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig8','b.png'));

% Plot GC thickness vs axial length, ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick))]
h=regressionPlot(comboTable.Axial_Length_average, gcMeanThick_reconstruct', 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig8','c.png'));

% Plot GC thickness vs axial length, ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick))]
h=regressionPlot(comboTable.Axial_Length_average, gcMeanThick_reconstruct_orig', 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig8','d.png'));