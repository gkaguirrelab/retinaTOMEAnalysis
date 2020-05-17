for ii = 1:50
    profileFit = GCVolPCAScoreExpandedSmoothed(:,1:nDimsToUse)*GCVolPCACoeff(ii,1:nDimsToUse)';
    profileFit(isnan(gcVolumePerDegSq(:,ii)))=nan;
    meanFitGCVol(ii) = nanmean(profileFit);
    
    profileFit = GCVolPCAScoreExpandedSmoothed(:,1:nDimsToUse)*adjustedGCVolPCACoeff(ii,1:nDimsToUse)';
    profileFit(isnan(gcVolumePerDegSq(:,ii)))=nan;
    meanAdjustedGCVol(ii) = nanmean(profileFit);

end

% Get the mean fit volumes and add these to the combo table
fitVolTable = cell2table([num2cell(str2double(subList)'),num2cell(meanFitGCVol)',num2cell(meanAdjustedGCVol)'],...
    'VariableNames',{'AOSO_ID','meanFitGCVol','meanAdjustedGCVol'});
comboTable = join(comboTable,fitVolTable,'Keys','AOSO_ID');

% Load the post-retinal anatomy table and make a joined table
opts = detectImportOptions(p.Results.anatMeasuresFileName);
anatMeasuresTable = readtable(p.Results.anatMeasuresFileName, opts);
anatComparisonTable = join(anatMeasuresTable,comboTable,'Keys','AOSO_ID');

%% Model optic chiasm volume
y = anatComparisonTable.Optic_Chiasm;

% Create an X model
X = [anatComparisonTable.meanAdjustedGCVol, anatComparisonTable.meanFitGCVol, anatComparisonTable.SupraTentorialVol./1e9];
mdl = fitlm(X,y,'linear');

X = [anatComparisonTable.meanAdjustedGCVol, anatComparisonTable.meanFitGCVol, anatComparisonTable.SupraTentorialVol./1e9, ones(size(anatComparisonTable.meanAdjustedGCVol))];
b = regress(y,X);
yFit = X*b;
