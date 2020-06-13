function fig10_postRetinalAnatomy(p,GCVolPCACoeff,gcVolumePerDegSq,adjustedGCVolPCACoeff,comboTable,GCVolPCAScoreExpandedSmoothed,nDimsToUse,subList,saveDir)


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


%% Report the correlation of optic chiasm volume with some measures
measureSet = {'gcMeanThick','meanFitGCVol','meanAdjustedGCVol','Height_inches','Weight_pounds','Age','Axial_Length_average','Gender'};
y = anatComparisonTable.Optic_Chiasm;
nBoots = 1000;
for ii = 1:length(measureSet)
    if strcmp(measureSet{ii},'Gender')
        tmp = anatComparisonTable.(measureSet{ii});
        x = zeros(size(tmp));
        x(strcmp(tmp,'M'))=1;        
    else
        x = anatComparisonTable.(measureSet{ii});
    end
    [R,P] = corrcoef(y,x);
    for bb = 1:nBoots
        bootSamp = randsample(length(y),length(y),true);
        bootR(bb) = corr(y(bootSamp),x(bootSamp));
    end
    str = sprintf(['Correlation of optic chiasm volume with ' measureSet{ii} ' = %2.2f ± %2.2f (sem), p = %2.5f \n'],R(1,2),std(bootR),P(1,2));
    fprintf(str);
end


%% Model optic chiasm volume

% Create an X model
X = [anatComparisonTable.meanAdjustedGCVol, anatComparisonTable.meanFitGCVol];
X(:,1) = X(:,1)-mean(X(:,1));
X(:,2) = X(:,2)-mean(X(:,2));
X(:,1) = X(:,1)/2;
X(:,2) = X(:,2)/2;

mdl = fitlm(X,y,'linear');

figHandle = figure();
h = mdl.plot;
h(1).Marker = 'o';
h(1).MarkerEdgeColor = 'none';
h(1).MarkerFaceColor = 'r';

h(2).Color = [0.5 0.5 0.5];
h(3).Color = [0.5 0.5 0.5];

xlabel('Modeled deviation from mean GC Tissue Volume [mm^3 / deg^2]')
ylabel('Optic chiasm volume [mm^3]')

setTightFig
saveas(figHandle,fullfile(saveDir,'fig10','b.pdf'));

end