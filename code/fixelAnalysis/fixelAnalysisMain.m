function fixelAnalysisMain(varargin)
% Relate variation in GCL volue to fixel measurements on the visual pathway
%
% Examples:
%{
%}

%% Set the dropboxBaseDir
% We need this for the default loations of some the directories
dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));
%dropboxBaseDir ='C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)';


%% Parse vargin
p = inputParser;

% Optional analysis params
p.addParameter('showPlots',true,@islogical);
p.addParameter('horizVolFile',fullfile(dropboxBaseDir,'AOSO_analysis','GCPaperFigures_horiz','gcVolumeData.mat'),@ischar);
p.addParameter('vertVolFile',fullfile(dropboxBaseDir,'AOSO_analysis','GCPaperFigures_vert','gcVolumeData.mat'),@ischar);
p.addParameter('eyeModelsFileName',fullfile(dropboxBaseDir,'AOSO_analysis','eyeModels','eyeModels.mat'),@ischar);
p.addParameter('anatMeasuresFileName',fullfile(getpref('retinaTOMEAnalysis','projectBaseDir'),'data','visualPathwayAnatMeasures.xlsx'),@ischar);
p.addParameter('fixelDataDir',fullfile(getpref('retinaTOMEAnalysis','projectBaseDir'),'data','fixelResults'),@ischar);
p.addParameter('smoothPCAFactor',0.85,@isscalar);

% Check the parameters
p.parse(varargin{:});

% Load the horizontal and vertical gc tissue volume
load(p.Results.vertVolFile,'gcVolumePerDegSq','badIdx','XPos_Degs','subList');
gcVolumePerDegSq_vert = gcVolumePerDegSq;
badIdx_vert = badIdx;
XPos_Degs_vert = XPos_Degs;
subList_vert = subList;

load(p.Results.horizVolFile,'gcVolumePerDegSq','badIdx','XPos_Degs','comboTable','subList');
gcVolumePerDegSq_horiz = gcVolumePerDegSq;
badIdx_horiz = badIdx;
XPos_Degs_horiz = XPos_Degs;
subList_horiz = subList;

if isequal(subList_vert,subList_horiz)
    subList = subList_vert;
else
    error('These should match')
end

nDimsToUse = 6; % number of PCA components to use.


% Concatenate the horizontal and vertical for the purposes of the PCA
% analysis and axial length correction
gcVolumePerDegSq = [gcVolumePerDegSq_horiz; gcVolumePerDegSq_vert];
badIx = [badIdx_horiz;badIdx_vert];
posShift = 2*max(XPos_Degs_horiz)+1;
XPos_Degs = [XPos_Degs_horiz; XPos_Degs_vert+posShift] ;


%% Conduct PCA upon the tissue volume data
% Perform the PCA and smooth components
[GCVolPCAScoreExpanded, GCVolPCAScoreExpandedSmoothed, GCVolPCACoeff, GCVolPCAVarExplained] = createVolumePCA(gcVolumePerDegSq,badIdx,XPos_Degs,nDimsToUse,p.Results.smoothPCAFactor,'both');

% Adjust each coeff with the axial length contribution, also create some
% synthetic ones
[adjustedGCVolPCACoeff,synGCVolPCACoeff] = adjustAndSynthPCAWithAxialLength(nDimsToUse,GCVolPCACoeff,comboTable.Axial_Length_average);

% Get the mean fit and adjusted GC volume
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


% Instantiate a flywheel object
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download the fixel results
% To find these, get the ID for the session (which is in the URL of the web
% GUI, and then use this command to get a list of the analyses associated
% with that session, and then find the analysis ID the we want.
%
%{
    projectName = 'flywheelMRSupport';
    fw = flywheel.Flywheel(getpref(projectName,'flywheelAPIKey'));
    sessionID = '6036e964070e35d8850b67c9';
    analysisList = fw.getSessionAnalyses(sessionID);
%}
% Right then left optic tract
laterality = {'right','left'};
analysisIDs = {'60525004a51b90af39c85125','60524fe6620334967cd8daf3'};
fileNames = {'fc_stats.csv','fd_stats.csv','fdc_stats.csv'};
for ll = 1:length(laterality)
    for ff = 1:length(fileNames)
        saveName = fullfile(p.Results.fixelDataDir,[laterality{ll} '_' fileNames{ff}]);
        fw.downloadOutputFromAnalysis(analysisIDs{ll},fileNames{ff},saveName);
        
        % Now load the file
        opts = detectImportOptions(saveName);
        fixelData = readtable(saveName, opts);
        if ll==1 && ff==1
            fixelTable = fixelData(:,1:2);
            fixelTable.Properties.VariableNames{2} = [laterality{ll} '_' fileNames{ff}(1:3)];
        else
            subTable = fixelData(:,1:2);
            subTable.Properties.VariableNames{2} = [laterality{ll} '_' fileNames{ff}(1:3)];
            fixelTable=join(fixelTable,subTable);
        end
    end
end

% Massage the fixelTable to match up with the comboTable
fixelTable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable.TOME_ID = strrep(fixelTable.TOME_ID,'fod_','');
fixelSet = {'fc_','fd_','fdc'};
for ff = 1:length(fixelSet)
    % Report the correlation of left and right
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    [R,P] = corrcoef(fixelValR,fixelValL);
    fprintf([fixelSet{ff} ' correlation left with right: %2.2f \n'],R(1,2));
    fixelTable.(fixelSet{ff}) = mean([fixelValR, fixelValL],2);
end

fixelComparisonTable = join(comboTable(ismember(comboTable.TOME_ID,fixelTable.TOME_ID),:),fixelTable,'Keys','TOME_ID');


%% Report the correlation of fixel values with RGC values
measureSet = {'gcMeanThick','meanFitGCVol','meanAdjustedGCVol','Height_inches','Weight_pounds','Age','Axial_Length_average','Gender'};

for ff = 1:length(fixelSet)
    y = fixelComparisonTable.(fixelSet{ff});
    nBoots = 1000;
    for ii = 1:length(measureSet)
        if strcmp(measureSet{ii},'Gender')
            tmp = fixelComparisonTable.(measureSet{ii});
            x = zeros(size(tmp));
            x(strcmp(tmp,'M'))=1;
        else
            x = fixelComparisonTable.(measureSet{ii});
        end
        [R,P] = corrcoef(y,x);
        for bb = 1:nBoots
            bootSamp = randsample(length(y),length(y),true);
            bootR(bb) = corr(y(bootSamp),x(bootSamp));
        end
        str = sprintf(['Correlation of ' fixelSet{ff} ' with ' measureSet{ii} ' = %2.2f ± %2.2f (sem), p = %2.5f \n'],R(1,2),std(bootR),P(1,2));
        fprintf(str);
    end
end    


%% Model fc by GC values

y = log10(fixelComparisonTable.fc_);

% Create an X model
X = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanFitGCVol];
X(:,1) = X(:,1)-mean(X(:,1));
X(:,2) = X(:,2)-mean(X(:,2));
X(:,1) = X(:,1)/2;
X(:,2) = X(:,2)/2;

mdl = fitlm(X,y,'linear')

figHandle = figure();
h = mdl.plot;
h(1).Marker = 'o';
h(1).MarkerEdgeColor = 'none';
h(1).MarkerFaceColor = 'r';

h(2).Color = [0.5 0.5 0.5];
h(3).Color = [0.5 0.5 0.5];

xlabel('Modeled deviation from mean GC Tissue Volume [mm^3 / deg^2]')
ylabel('optic tract fc')

setTightFig

end

