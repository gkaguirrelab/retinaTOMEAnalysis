function fig9_SynthProfileAndALRelationship(XPos_Degs,scoreExpandedSmoothed,synCoeff,ALRange,nDimsToUse,saveDir)
% Plot the synthesized reconstructions by axial length
gcVec_reconstruct = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(:,1:nDimsToUse)';
meangcVec_reconstruct = mean(gcVec_reconstruct,2);
gcMedianThick_reconstruct = nanmedian(gcVec_reconstruct,1);
gcMeanThick_reconstruct = nanmean(gcVec_reconstruct,1);


% Plot the GC thickness functions, ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))]
h=profilePlot(XPos_Degs, gcVec_reconstruct, meangcVec_reconstruct, 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],1);
setTightFig
saveas(h,fullfile(saveDir,'fig9','a.png'));

