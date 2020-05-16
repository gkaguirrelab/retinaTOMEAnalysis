function fig4_smoothedPCAReconstruction(gcVolumePerDegSq,scoreExpandedSmoothed,coeff,nDimsToUse,saveDir)
h=figure;
set(gcf,'color','w');
for ii = 1:49
    subplot(7,7,ii);
    plot(gcVolumePerDegSq(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff(ii,1:nDimsToUse)';
    profileFit(isnan(gcVolumePerDegSq(:,ii)))=nan;
    plot(profileFit,'-r','LineWidth',1);
    axis off
end
suptitle('Original and Reconstructed GC Tissue Volume Profiles by Subject')
setTightFig
saveas(h,fullfile(saveDir,'fig4','a.pdf'));
