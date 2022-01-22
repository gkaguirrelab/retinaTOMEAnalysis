
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsMm.mat');
load(individualFitFile,'RSquaredSet');
mmRSquaredSet = RSquaredSet;

sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsDeg.mat');
load(individualFitFile,'RSquaredSet');
degRSquaredSet = RSquaredSet;

[~,p,~,stats] = ttest(mmRSquaredSet(4,:),degRSquaredSet(4,:));

[p,~,stats] = ranksum(mmRSquaredSet(2,:),degRSquaredSet(2,:));
