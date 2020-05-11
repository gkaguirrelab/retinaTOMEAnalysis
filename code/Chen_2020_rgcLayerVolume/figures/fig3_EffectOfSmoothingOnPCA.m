function fig3_EffectOfSmoothingOnPCA(explained,scoreExpanded,scoreExpandedSmoothed,nDimsToUse,saveDir)
% Show the effect of smoothing on the PCA scores
h=figure;
set(gcf,'color','w');
shapeExplained = explained ./ sum(explained(2:end));
shapeExplained(1) = nan;
for ii = 1:nDimsToUse
    subplot(3,2,ii);
    plot(scoreExpanded(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    plot(scoreExpandedSmoothed(:,ii),'-r','LineWidth',1);
    str = sprintf('PC%d, shape var explained: %2.2f',ii,shapeExplained(ii));
    title(str);
    axis off
end
suptitle('Original and Smoothed PCA Components')
setTightFig
saveas(h,fullfile(saveDir,'fig3','a.png'));
