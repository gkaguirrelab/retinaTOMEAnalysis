% Plot the reconstructions with the adjustment
h=figure
set(gcf,'color','w');
[ALsorted, ALsortedIndx] = sort(axialLengths);
counter=0;
for ii = ALsortedIndx
    counter = counter+1;
    subplot(8,7,counter);
    plot(gcVolumePerDegSq(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*adjustedCoeff(ii,1:nDimsToUse)';
    plot(profileFit,'-r','LineWidth',1);
    ylim([-1 8])
    axis off
end
setTightFig
saveas(h,fullfile(saveDir,'fig6','a.png'));

%suptitle('Original and AL influence adjusted gc tissue volume profiles by subject')

