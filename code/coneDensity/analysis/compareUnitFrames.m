% We can consider the model in terms of mm of retina, or degrees of visual
% angle. The model was fit to each subject in each of the two unit frames.
% We can test if the R2 of the model fit is consistently higher in one
% frame or the other.

% Load the R2 results for the deg unit frame
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsDeg.mat');
load(individualFitFile,'RSquaredSet');
degRSquaredSet = RSquaredSet;

% Load the R2 results for the mm unit frame (and some other files for
% later)
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsMm.mat');
load(individualFitFile,'RSquaredSet');
mmRSquaredSet = RSquaredSet;

% Conduct a Wilcoxon rank-sum test that asks if the mm frame is better than
% the degree frame for just the eccentricity model
[p,~,stats] = ranksum(mmRSquaredSet(2,:),degRSquaredSet(2,:));
fprintf('Test for higher R2 values for eccentricity model in the mm as compared to the deg framework: Wilcoxon z= %2.2f, p= %2.2f \n',stats.zval,p)

[p,~,stats] = ranksum(mmRSquaredSet(4,:),degRSquaredSet(4,:));
fprintf('Test for higher R2 values for polar angle model in the mm as compared to the deg framework: Wilcoxon z= %2.2f, p= %2.4f \n',stats.zval,p)

% The result is that the double-exponential model fits the data about as
% well when cast in mm or degrees. However, the variation in polar angle is
% better modeled in mm. This suggests that different eyes have a similar
% polar angle variation when cast in mm, perhaps because this variation
% across the retina is locked to physical distance from the fovea, and does
% not scale proportinately in eyes of different sizes.

figure
rSquareTitles = {'RSquaredFull','RSquaredExponentialOnly','RSquaredPolarOnly','RSquaredPolarResiduals'};
indexSet = [2 4];
for ii=1:length(indexSet)
    subplot(1,2,ii)
    [~,edges] = histcounts(mmRSquaredSet(indexSet(ii),:));
    binsize = edges(2)-edges(1);
    edges = 0:binsize:1+(binsize/2);
    centers = binsize/2:binsize:1;
    centerI = 0:0.01:1;
    N = histcounts(mmRSquaredSet(indexSet(ii),:),edges);
    plot(centers,N,'ok'); hold on
    Nq = interp1(centers,N,centerI,'makima',0);
    plot(centerI,Nq,'-r');
    xlabel('R squared')
    ylabel('counts [subjects]')
    title(rSquareTitles{indexSet(ii)})
end


