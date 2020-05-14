function fig9_SynthProfileAndALRelationship(XPos_Degs,comboTable,scoreExpandedSmoothed,synCoeff,nDimsToUse,saveDir)
% Plot the synthesized reconstructions by axial length
[~, myopeIdx] = min(comboTable.Axial_Length_average);
[~, hyperopeIdx] = max(comboTable.Axial_Length_average);
[~, emmetropeIdx] = min( abs(comboTable.Axial_Length_average - 23.58) );
idx = [myopeIdx, hyperopeIdx, emmetropeIdx];
gcVec_reconstruct = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(idx,1:nDimsToUse)';

% Plot the GC volume functions
str = sprintf('Synthesized GC vol profiles for AL = %2.2f, %2.2f, %2.2g,',comboTable.Axial_Length_average(idx));
h=profilePlot(XPos_Degs, gcVec_reconstruct, [], 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig9','a.png'));

% % Convert from volume to thickness
% mmSqPerDegSq(:,ss) = mmPerDegPolyFit{idx}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
%     gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);
% 
% mmPerDegPolyFit{idx}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
%     gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);