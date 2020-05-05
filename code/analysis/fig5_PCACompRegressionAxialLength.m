for i = 1:nDimsToUse
h=regressionPlot(comboTable.Axial_Length_average, coeff(:,i), 'Axial Length [mm]',['PC' num2str(i) ' Loading'], [],p.Results.showPlots);
setTightFig
saveas(h,fullfile(saveDir,'fig5',[num2str(i) '.png']));
end