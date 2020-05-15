function fig2_ThickAndVolRelationships(XPos_Degs, gcVec, meanGCVec, mmSqPerDegSq, gcVolumePerDegSq, meanGCVolumePerDegSqProfile,comboTable,saveDir)
% Plot the GC thickness functions, ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, gcVec, meanGCVec, 'Eccentricity [deg visual angle]','GC Thickness [mm]',[],1);
setTightFig
saveas(h,fullfile(saveDir,'fig2','a.png'));

% Plot GC thickness vs axial length
str = ['Fig 2b: Axial length vs. mean GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick)) '\n'];
fprintf(str);
h=regressionPlot(comboTable.Axial_Length_average, comboTable.gcMeanThick, 'Axial Length [mm]','Mean GC Thickness [mm]', [],1);
setTightFig
saveas(h,fullfile(saveDir,'fig2','b.png'));

% Plot the mmSqPerDegSq functions
h=profilePlot(XPos_Degs, mmSqPerDegSq, nanmean(mmSqPerDegSq,2), 'Eccentricity [deg visual angle]','[mm^2 / deg^2]', [],1);
setTightFig
saveas(h,fullfile(saveDir,'fig2','d.png'));

% Plot gc tissue volume profile
h=profilePlot(XPos_Degs, gcVolumePerDegSq, meanGCVolumePerDegSqProfile, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]', [] ,1);
setTightFig
saveas(h,fullfile(saveDir,'fig2','e.png'));

% Plot mean GC tissue volume vs axial length
str = ['Fig 2f: Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq)) '\n'];
fprintf(str);
h=regressionPlot(comboTable.Axial_Length_average, comboTable.gcVolumePerDegSq, 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', ...
    [],1);
setTightFig
saveas(h,fullfile(saveDir,'fig2','f.png'));