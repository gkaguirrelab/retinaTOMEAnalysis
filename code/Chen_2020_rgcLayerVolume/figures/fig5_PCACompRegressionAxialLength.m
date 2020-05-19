function fig5_PCACompRegressionAxialLength(axialLengths,coeff,nDimsToUse,saveDir)

for i = 1:nDimsToUse
h=regressionPlot(axialLengths, coeff(:,i), 'Axial Length [mm]',['PC' num2str(i) ' Loading'], [],1);
setTightFig
saveas(h,fullfile(saveDir,'fig5',[num2str(i) '.pdf']));
end