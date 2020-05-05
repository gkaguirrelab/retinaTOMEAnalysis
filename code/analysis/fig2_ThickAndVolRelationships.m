% Plot the GC thickness functions, ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, gcVec, meanGCVec, 'Eccentricity [deg visual angle]','GC Thickness [microns]',[],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig2','a.png'));

% Plot GC thickness vs axial length, ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick))]
h=regressionPlot(comboTable.Axial_Length_average, comboTable.gcMeanThick, 'Axial Length [mm]','Mean GC Thickness [microns]', [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig2','b.png'));

% Plot the mmSqPerDegSq functions, ['mm square retina per degree square visuall angle (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, mmSqPerDegSq, nanmean(mmSqPerDegSq,2), 'Eccentricity [deg visual angle]','[mm^2 / deg^2]', [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig2','d.png'));

% Plot gc tissue volume profile, ['GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, gcVolumePerDegSq, meanGCVolumePerDegSqProfile, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]', [] ,p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig2','r.png'));

% Plot mean GC tissue volume vs axial length, ['Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq))]
h=regressionPlot(comboTable.Axial_Length_average, comboTable.gcVolumePerDegSq, 'Axial Length [mm]','Mean GC Tissue Volume [mm^3 / deg^2]', ...
    [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig2','f.png'));