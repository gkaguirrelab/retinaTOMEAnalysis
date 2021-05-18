function fixelAnalysisMain(varargin)
% Relate variation in GCL volue to fixel measurements on the visual pathway
%
% Examples:
%{
%}

%% Set the dropboxBaseDir
% We need this for the default loations of some the directories
dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));
% dropboxBaseDir='C:\Users\ozenc\Dropbox (Aguirre-Brainard Lab)';

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
% fw = flywheel.Flywheel('upenn.flywheel.io:DTIiZcuXBVlpJmCLZt');

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
analysisIDs = {'60a33c76b4a131197e7bfaa8','60a33c7617fcfbb03ffeacf6'};
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

% Sort rows by subject ID, so that it will be easier to add other measures
fixelTable = sortrows(fixelTable);

% Add FA and MD to the fixeltable 
laterality = {'right','left','right','left'};
analysisIDs = {'60a3517044aedd9b5ef5e603','60a3514a50e3c45f1bfead13','60a3510e038ae96d18f5e69c','60a350c8c711ebba1ca23c44'};
fileNames = {'FA_stats.csv', 'FA_stats.csv', 'MD_stats.csv', 'MD_stats.csv'};
for ll = 1:length(laterality)
    saveName = fullfile(p.Results.fixelDataDir,[laterality{ll} '_' fileNames{ll}]);
    fw.downloadOutputFromAnalysis(analysisIDs{ll},fileNames{ll},saveName);
        
    % Now load the file
    opts = detectImportOptions(saveName);
    dtiData = readtable(saveName, opts);
        
    % Sort the table
    dtiData = sortrows(dtiData);
        
    % Add 
    subTable = dtiData(:,1:2);
    subTable.Properties.VariableNames{1} = 'TOME_ID';
    subTable.Properties.VariableNames{2} = [laterality{ll} '_' fileNames{ll}(1:2)];
    if strcmp(fixelTable{1,1}{1}(1:7), 'preproc')
        for ii = 1:height(fixelTable)
            fixelTable{ii,1}{1} = fixelTable{ii,1}{1}(12:end);
        end
    end
    fixelTable=join(fixelTable,subTable);
end

%% Collect total intracranial volume from Freesurfer files 

% Get tome subjects
projects = fw.projects();
tome = projects{1,1};
subjects = tome.subjects();
subjectLength = length(subjects);

% Create a folder in fixelData dir for extracting aseg files from
% Freesurfer directories
asegDir = fullfile(p.Results.fixelDataDir, 'asegDir');
if ~exist(asegDir, 'dir')
    system(['mkdir' ' ' asegDir]);
end

% Create empty matrices for subject and intracranialVol
subjectNames = [];
intracranialVol = [];

% Loop through subjects
for sub = 1:subjectLength
    % If label starts with TOME and it's not TOME_3027 (because that
    % subject was discarted), process the subject
    if strcmp(subjects{sub}.label(1:4), 'TOME') && ~strcmp(subjects{sub}.label(6:end), '3027')      
        % Get subject name and save name where aseg files will be saved
        subject = subjects{sub,1};
        subjectLabel = subject.label;
        subjectNames = [subjectNames; {subjectLabel}];
        saveName = fullfile(asegDir,[subjectLabel '_aseg.stats']); 
        
        % If aseg file does not exist, download it.
        if ~isfile(saveName)
            % Get session
            sessions = subject.sessions();
            for ses = 1:length(sessions)
                session = sessions{ses,1};
                % Get analysis and loop through
                analyses = session.analyses();
                for a = 1:length(analyses)
                    % Find analyses thay contain freesurfer in name
                    if contains(analyses{a,1}.label, 'freesurfer')
                        freesurferAnalysisContainer = analyses{a,1};
                        analysisIdTag = freesurferAnalysisContainer.id;
                        zipFile = ['freesurfer-recon-all_' subjectLabel '_' analysisIdTag '.zip'];  
                        % Download only the aseg files from the whole zip
                        freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/stats/aseg.stats'], saveName);
                    end
                end
            end
        end
        % Load the aseg files for each subject, extract the intrcranial
        % volume and save it to the intracranial matrix 
        asegFileLoaded = textread(saveName, '%s');
        intraCranialVolume = str2num(asegFileLoaded{297});
        intracranialVol = [intracranialVol; intraCranialVolume];
    end
end

% Combine subject and intracranial volume in a table and sort rows.
intracranialTable = table(subjectNames, intracranialVol);
intracranialTable = sortrows(intracranialTable);
intracranialTable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,intracranialTable);

fixelSet = {'fc_','fd_','fdc', 'FA', 'MD'};
for ff = 1:length(fixelSet)
    % Report the correlation of left and right
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    [R,P] = corrcoef(fixelValR,fixelValL);
    fprintf([fixelSet{ff} ' correlation left with right: %2.2f \n'],R(1,2));
    fixelTable.(fixelSet{ff}) = mean([fixelValR, fixelValL],2);
end

% Adding intracranial volume to the fixel set here, so that it will also be
% correlated with measure set variables 
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'intracranialVol'};

fixelComparisonTable = join(comboTable(ismember(comboTable.TOME_ID,fixelTable.TOME_ID),:),fixelTable,'Keys','TOME_ID');


%% Report the correlation of fixel values with RGC values
measureSet = {'gcMeanThick','meanFitGCVol','meanAdjustedGCVol','Height_inches','Weight_pounds','Age','Axial_Length_average','Gender','intracranialVol'};

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

% Partial Correlations
fprintf('\nPartial correlations\n')
controlFor = {'Height_inches','Weight_pounds','intracranialVol'};
for ii = 1:length(controlFor)
    [rho, p] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_, fixelComparisonTable.(controlFor{ii}));
    fprintf(['Partial correlation FC with meanAdjustedGCVol controlled for ' controlFor{ii} ': ' 'rho:' num2str(rho) ', p:' num2str(p) '\n'])
end
for ii = 1:length(controlFor)
    [rho, p] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fd_, fixelComparisonTable.(controlFor{ii}));
    fprintf(['Partial correlation FD with meanAdjustedGCVol controlled for ' controlFor{ii} ': ' 'rho:' num2str(rho) ', p:' num2str(p) '\n']) 
end

genderMatrix = fixelComparisonTable.Gender;
for ii = 1:length(genderMatrix)
    if strcmp(genderMatrix{ii}, 'M')
        genderMatrix{ii} = 1;
    else
        genderMatrix{ii} = 0;
    end
end
        
sizeMatrixFC = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.intracranialVol cell2mat(genderMatrix)];
[rho, p] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_, sizeMatrixFC);
fprintf(['\nPartial correlation FC with meanAdjustedGCVol controlled for height,weight,ICV,gender: ' 'rho:' num2str(rho) ', p:' num2str(p) '\n'])

sizeMatrixFD = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.intracranialVol cell2mat(genderMatrix)];
[rho, p] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fd_, sizeMatrixFD);
fprintf(['\nPartial correlation FD with meanAdjustedGCVol controlled for height,weight,ICV,gender: ' 'rho:' num2str(rho) ', p:' num2str(p) '\n'])

% Make controlled correlation plot for FC
tableFCvsSize = fitlm(sizeMatrixFC, fixelComparisonTable.fc_);
FCResiduals = tableFCvsSize.Residuals.Pearson;

tableGCvsSize = fitlm(sizeMatrixFC, fixelComparisonTable.meanAdjustedGCVol);
GCResidualsFC = tableGCvsSize.Residuals.Pearson;

[R,p] = corrcoef(GCResidualsFC, FCResiduals);
figure;
% add first plot in 2 x 1 grid    
scatter(GCResidualsFC, FCResiduals, 'MarkerFaceColor', 'k');
xlabel ('GCResiduals');
ylabel('FCResiduals');
title('FC vs MeanGC controlled for height/weight/ICV/Gender')
box 'on'
axis square;
set(gca,'Ticklength',[0 0])
%white background
set(gcf,'color','w');
refline
theStringR = sprintf(['R=' ' ' num2str(R(1,2))], GCResidualsFC,  FCResiduals);
theStringP = sprintf(['P=' ' ' num2str(p(1,2))], GCResidualsFC,  FCResiduals);
text(1.5, -1.5, theStringR, 'FontSize', 10);
text(1.5, -1.8, theStringP, 'FontSize', 10);

% Make controlled correlation map for FD
tableFDvsSize = fitlm(sizeMatrixFD, fixelComparisonTable.fd_);
FDResiduals = tableFDvsSize.Residuals.Pearson;

tableGCvsSize = fitlm(sizeMatrixFD, fixelComparisonTable.meanAdjustedGCVol);
GCResidualsFD = tableGCvsSize.Residuals.Pearson;

[R,p] = corrcoef(GCResidualsFD, FDResiduals);
figure;
% add first plot in 2 x 1 grid    
scatter(GCResidualsFD, FDResiduals, 'MarkerFaceColor', 'k');
xlabel ('GCResiduals');
ylabel('FDResiduals');
title('FD vs MeanGC controlled for height/weight/ICV/Gender')
box 'on'
axis square;
set(gca,'Ticklength',[0 0])
%white background
set(gcf,'color','w');
refline
theStringR = sprintf(['R=' ' ' num2str(R(1,2))], GCResidualsFD,  FDResiduals);
theStringP = sprintf(['P=' ' ' num2str(p(1,2))], GCResidualsFD,  FDResiduals);
text(1.5, -2.5, theStringR, 'FontSize', 10);
text(1.5, -2.8, theStringP, 'FontSize', 10);
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

