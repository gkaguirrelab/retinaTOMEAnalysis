function fixelAnalysisMain(varargin)
% Relate variation in GCL volue to fixel measurements on the visual pathway
%
% Examples:
%{
%}

%% Set the dropboxBaseDir
% We need this for the default loations of some the directories
% dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));
dropboxBaseDir='/home/ozzy/Dropbox (Aguirre-Brainard Lab)';

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
% analysisIDs = {'613222eae4575d3ae7e439e1','613222e90cb137b99763de1d'};
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

% Right then left optic radiation
laterality = {'right','left'};
% analysisIDs = {'613dd3a1af9c5aae928f7b26','613dd3a150fe1777a6965bad'};
analysisIDs = {'60ed117b6d2438c15c96c6bd','60ed1153dcf573726496c77d'}; % NOT A SINGLE SHELL
fileNames = {'fc_stats.csv','fd_stats.csv','fdc_stats.csv'};
for ll = 1:length(laterality)
    for ff = 1:length(fileNames)
        opticRadiationFolder = fullfile(p.Results.fixelDataDir, 'opticRadiation');
        if ~exist(opticRadiationFolder, 'dir')
            system(['mkdir' ' ' opticRadiationFolder]);
        end
        saveName = fullfile(opticRadiationFolder,[laterality{ll} '_opticRadiation_' fileNames{ff}]);
        fw.downloadOutputFromAnalysis(analysisIDs{ll},fileNames{ff},saveName);
        
        % Now load the file
        opts = detectImportOptions(saveName);
        fixelData = readtable(saveName, opts);
        if ll==1 && ff==1
            fixelTableRadiation = fixelData(:,1:2);
            fixelTableRadiation.Properties.VariableNames{2} = [laterality{ll} '_' fileNames{ff}(1:3) 'opticRadiation'];
        else
            subTable = fixelData(:,1:2);
            subTable.Properties.VariableNames{2} = [laterality{ll} '_' fileNames{ff}(1:3) 'opticRadiation'];
            fixelTableRadiation=join(fixelTableRadiation,subTable);
        end
    end
end

% Massage the fixelTable to match up with the comboTable
fixelTableRadiation.Properties.VariableNames{1} = 'TOME_ID';
fixelTableRadiation.TOME_ID = strrep(fixelTableRadiation.TOME_ID,'fod_','');

% Sort rows by subject ID, so that it will be easier to add other measures
fixelTableRadiation = sortrows(fixelTableRadiation);
fixelTable = join(fixelTable, fixelTableRadiation);

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
subjectDataFolder = fullfile(p.Results.fixelDataDir, 'subjectDataFolder');
if ~exist(subjectDataFolder, 'dir')
    system(['mkdir' ' ' subjectDataFolder]);
end

% Create empty matrices for subject and intracranialVol
subjectNames = [];
intracranialVol = [];
leftLGN = [];
rightLGN = [];
new_leftLGN = [];
new_rightLGN = [];
V1surfaceLeft = [];
V1surfaceRight = [];
V1ThicknessLeft = [];
V1ThicknessRight = [];
opticChiasmVolume = [];
opticChiasmVolumeHCP = [];
intracranialHCP = [];

% Loop through subjects
for sub = 1:subjectLength
    % If label starts with TOME and it's not TOME_3027 (because that
    % subject was discarted), process the subject
    if strcmp(subjects{sub}.label(1:4), 'TOME') && ~strcmp(subjects{sub}.label(6:end), '3027')      
        % Get subject name and save name where aseg files will be saved
        subject = subjects{sub,1};
        subjectLabel = subject.label;
        subjectFolder = fullfile(subjectDataFolder, subjectLabel);
        if ~exist(subjectFolder, 'dir')
            system(['mkdir' ' ' subjectFolder]);
        end
        subjectNames = [subjectNames; {subjectLabel}];
        asegSaveName = fullfile(subjectFolder,[subjectLabel '_aseg.stats']);
        asegSaveNameHCP = fullfile(subjectFolder,[subjectLabel '_HCPaseg.stats']);
        lhWhite = fullfile(subjectFolder,[subjectLabel '_lh.white']);
        rhWhite = fullfile(subjectFolder,[subjectLabel '_rh.white']);        
        lhThickness = fullfile(subjectFolder,[subjectLabel '_lh.thickness']);
        rhThickness = fullfile(subjectFolder,[subjectLabel '_rh.thickness']);
        LGNSaveName = fullfile(subjectFolder,[subjectLabel '_ThalamicNuclei.v12.T1.volumes.txt']);
        new_LGNSaveName = fullfile(subjectFolder,[subjectLabel '_new_ThalamicNuclei.v12.T1.volumes.txt']);
        bayesPrfSaveNameLeft = fullfile(subjectFolder,['lh.' subjectLabel '_inferred_varea.mgz']);
        bayesPrfSaveNameRight = fullfile(subjectFolder,['rh.' subjectLabel '_inferred_varea.mgz']);
        
        % Do the next block if any of the stat files do not exist in the path
        if ~isfile(asegSaveName) || ~isfile(LGNSaveName) || ~isfile(new_LGNSaveName) || ~isfile(bayesPrfSaveNameLeft) || ~isfile(lhWhite) || ~isfile(lhThickness)
            sessions = subject.sessions();
            for ses = 1:length(sessions)
                session = sessions{ses,1};
                % Get analysis and loop through
                analyses = session.analyses();
                for a = 1:length(analyses)
                    % Get aseg from latest freesurfer runs     
                    if ~isfile(asegSaveName)
                        if contains(analyses{a,1}.label, 'freesurfer')
                            freesurferAnalysisContainer = analyses{a,1};
                            analysisTag = freesurferAnalysisContainer.id;
                            zipFile = ['freesurfer-recon-all_' subjectLabel '_' analysisTag '.zip'];  
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/stats/aseg.stats'], asegSaveName);
                        end
                    end
                        
                    % Find analyses that contain hcp-struct in the name to
                    % get thickness and white surfaces
                    if ~isfile(lhWhite) || ~isfile(rhWhite) || ~isfile(lhThickness)
                        if contains(analyses{a,1}.label, 'hcp-struct')  
                            freesurferAnalysisContainer = analyses{a,1};
                            zipFile = [subjectLabel '_hcpstruct.zip'];  
                            % Download only the aseg files from the whole zip                          
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/T1w/' subjectLabel  '/surf/lh.white'], lhWhite);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/T1w/' subjectLabel  '/surf/rh.white'], rhWhite);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/T1w/' subjectLabel  '/surf/lh.thickness'], lhThickness);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/T1w/' subjectLabel  '/surf/rh.thickness'], rhThickness);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, [subjectLabel '/T1w/' subjectLabel  '/stats/aseg.stats'], asegSaveNameHCP);                            
                        end
                    end
                    % Find analyses that contain segmentThalamicNuclei in
                    % the name
                    if ~isfile(LGNSaveName)
                        if contains(analyses{a,1}.label, 'segmentThalamicNuclei')
                            freesurferAnalysisContainer = analyses{a,1};
                            zipFile = [subjectLabel '.zip'];
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, 'ThalamicNuclei.v12.T1.volumes.txt', LGNSaveName);
                        end
                    end
                    % New LGN
                    if ~isfile(new_LGNSaveName)
                        if contains(analyses{a,1}.label, 'HCP - segment-thalamic-nuclei')
                            freesurferAnalysisContainer = analyses{a,1};
                            zipFile = [subjectLabel '.zip'];
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, 'ThalamicNuclei.v12.T1.volumes.txt', new_LGNSaveName);
                        end
                    end
                    % Find analyses that contain bayesPRF in the name
                    if ~isfile(bayesPrfSaveNameLeft)
                        if contains(analyses{a,1}.label, 'bayesprf')
                            freesurferAnalysisContainer = analyses{a,1};
                            zipFile = [subjectLabel '_inferred_surface.zip'];
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, ['lh.' subjectLabel '_inferred_varea.mgz'], bayesPrfSaveNameLeft);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, ['rh.' subjectLabel '_inferred_varea.mgz'], bayesPrfSaveNameRight);
                        end
                    end
                end
            end
        end

        % Load the aseg files for each subject, extract the intrcranial
        % volume and save it to the intracranial matrix 
        asegFileLoaded = textread(asegSaveName, '%s');
        intraCranialVolume = str2num(asegFileLoaded{297});
        intracranialVol = [intracranialVol; intraCranialVolume];
        
        % Load the aseg files from HCP for comparison
        asegHCPLoaded = textread(asegSaveNameHCP, '%s');
        intraCranialHCPRaw = str2num(asegHCPLoaded{309});
        intracranialHCP = [intracranialHCP; intraCranialHCPRaw]; 
        
        % extract the optic chiasm from FS and HCP
        opticCh = str2num(asegFileLoaded{917});
        opticChiasmVolume = [opticChiasmVolume; opticCh];
        opticChHCP = str2num(asegHCPLoaded{929});
        opticChiasmVolumeHCP = [opticChiasmVolumeHCP; opticChHCP];         
        
        % Load Thalamic segmentation text and extract LGN values
        fid = fopen(LGNSaveName);
        leftLGNLine = fgetl(fid);
        leftLGNVal = str2num(leftLGNLine(10:end));
        rightLGNLine = fgetl(fid);
        rightLGNVal = str2num(rightLGNLine(11:end));
        leftLGN = [leftLGN; leftLGNVal];
        rightLGN = [rightLGN; rightLGNVal];
        fclose(fid);
       
        % Load Thalamic segmentation text and extract LGN values
        new_fid = fopen(new_LGNSaveName);
        new_leftLGNLine = fgetl(new_fid);
        new_leftLGNVal = str2num(new_leftLGNLine(10:end));
        new_rightLGNLine = fgetl(new_fid);
        new_rightLGNVal = str2num(new_rightLGNLine(11:end));
        new_leftLGN = [new_leftLGN; new_leftLGNVal];
        new_rightLGN = [new_rightLGN; new_rightLGNVal];
        fclose(new_fid);
            
        % Left V1 surface area and thickness
        [vertLeft,faceLeft] = freesurfer_read_surf(lhWhite);
        [vareaMapLeft, M, mr_parms, volsz] = load_mgh(bayesPrfSaveNameLeft);
        vareaMapLeft = squeeze(vareaMapLeft);
        vertIdxLeft = find(vareaMapLeft==1);
        V1SurfaceAreaLeft = calcSurfaceArea(vertLeft,faceLeft,vertIdxLeft);     
        V1surfaceLeft = [V1surfaceLeft; V1SurfaceAreaLeft];
        thicknessLeft = read_curv(lhThickness);
        thicknessLeft = thicknessLeft(vertIdxLeft);
        V1ThicknessLeft = [V1ThicknessLeft; mean(thicknessLeft)];       
        
        % Right V1 surface area and thickness
        [vertRight,faceRight] = freesurfer_read_surf(rhWhite);
        [vareaMapRight, M, mr_parms, volsz] = load_mgh(bayesPrfSaveNameRight);         
        vareaMapRight = squeeze(vareaMapRight);
        vertIdxRight = find(vareaMapRight==1);
        V1SurfaceAreaRight = calcSurfaceArea(vertRight,faceRight,vertIdxRight);     
        V1surfaceRight = [V1surfaceRight; V1SurfaceAreaRight];
        thicknessRight = read_curv(rhThickness);
        thicknessRight = thicknessRight(vertIdxRight);
        V1ThicknessRight = [V1ThicknessRight; mean(thicknessRight)];  
    end
end

% Combine subject and intracranial volume in a table and sort rows.
intracranialTable = table(subjectNames, intracranialVol);
intracranialTable = sortrows(intracranialTable);
intracranialTable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,intracranialTable);

% HCP intracranial
intracranialTableHCP = table(subjectNames, intracranialHCP);
intracranialTableHCP = sortrows(intracranialTableHCP);
intracranialTableHCP.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,intracranialTableHCP);

% FS optic chiasm
chiasmTableFS = table(subjectNames, opticChiasmVolume);
chiasmTableFS = sortrows(chiasmTableFS);
chiasmTableFS.Properties.VariableNames{1} = 'TOME_ID';
chiasmTableFS.Properties.VariableNames{2} = 'opticChiasmVolume_FS';
fixelTable=join(fixelTable,chiasmTableFS);

% Combine subject and HCP optic chiasm volume in a table and sort rows.
chiasmTableHCP = table(subjectNames, opticChiasmVolumeHCP);
chiasmTableHCP = sortrows(chiasmTableHCP);
chiasmTableHCP.Properties.VariableNames{1} = 'TOME_ID';
chiasmTableHCP.Properties.VariableNames{2} = 'opticChiasmVolume_HCP';
fixelTable=join(fixelTable,chiasmTableHCP);

% % Get Min's optic chiasm results
% Optic_Chiasm_xls = readtable(fullfile(p.Results.fixelDataDir, 'MinvisualPathwayAnatMeasures.xlsx'));
% Optic_Chiasm = table(Optic_Chiasm_xls.TOME_ID, Optic_Chiasm_xls.Optic_Chiasm);
% Optic_Chiasm.Properties.VariableNames{1} = 'TOME_ID';
% Optic_Chiasm.Properties.VariableNames{2} = 'opticChiasmVolume_Min';
% fixelTable=join(fixelTable,Optic_Chiasm);

% Combine subject and LGN in a table and sort rows.
LGN = (leftLGN + rightLGN) / 2;
LGNTable = table(subjectNames, leftLGN, rightLGN, LGN);
LGNTable = sortrows(LGNTable);
LGNTable.Properties.VariableNames{1} = 'TOME_ID';
LGNTable.Properties.VariableNames{2} = 'left_LGN';
LGNTable.Properties.VariableNames{3} = 'right_LGN';
fixelTable=join(fixelTable,LGNTable);

% Combine subject and LGN in a table and sort rows.
new_LGN = (new_leftLGN + new_rightLGN) / 2;
new_LGNTable = table(subjectNames, new_leftLGN, new_rightLGN, new_LGN);
new_LGNTable = sortrows(new_LGNTable);
new_LGNTable.Properties.VariableNames{1} = 'TOME_ID';
new_LGNTable.Properties.VariableNames{2} = 'left_LGN_HCP';
new_LGNTable.Properties.VariableNames{3} = 'right_LGN_HCP';
fixelTable=join(fixelTable,new_LGNTable);

% Combine subject and V1 volume in a table and sort rows.
V1surface = (V1surfaceLeft + V1surfaceRight) / 2;
V1SurfaceTable = table(subjectNames, V1surface);
V1SurfaceTable = sortrows(V1SurfaceTable);
V1SurfaceTable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,V1SurfaceTable);

% Combine subject and V1 volume in a table and sort rows.
V1thickness = (V1ThicknessLeft + V1ThicknessRight) / 2;
V1ThicknessTable = table(subjectNames, V1thickness);
V1ThicknessTable = sortrows(V1ThicknessTable);
V1ThicknessTable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,V1ThicknessTable);

% Get SPM intracranial volume
SPMTIV = fullfile(p.Results.fixelDataDir, 'TIV.txt');
load(SPMTIV)
subSorted = sortrows(subjectNames);
TIVtable = table(subSorted, TIV);
TIVtable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,TIVtable);
%% Correlation of left right
fprintf('\n<strong>Correlation of left right hemisphere measurements of diffusion variables\n</strong>')
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation'};
for ff = 1:length(fixelSet)
    % Report the correlation of left and right
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    [R,P] = corrcoef(fixelValR,fixelValL);
    fprintf([fixelSet{ff} ' correlation left with right: %2.2f \n'],R(1,2));
    fixelTable.(fixelSet{ff}) = mean([fixelValR, fixelValL],2);
end

%% Report the correlation of fixel values with RGC values
% Variables to compare 
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation', 'intracranialVol', 'intracranialHCP', 'TIV', 'LGN', 'meanAdjustedGCVol', 'V1surface', 'V1thickness'};
fixelComparisonTable = join(comboTable(ismember(comboTable.TOME_ID,fixelTable.TOME_ID),:),fixelTable,'Keys','TOME_ID');

% Variables to compare against
measureSet = {'gcMeanThick','meanFitGCVol','meanAdjustedGCVol','Height_inches','Weight_pounds','Age','Axial_Length_average','Gender','intracranialVol', 'intracranialHCP', 'TIV', 'LGN', 'V1surface', 'V1thickness'};

% Names for size plots
names = {'Optic Tract FC','Optic Tract FD','fdc', 'FA', 'MD', 'Optic Radiation FC', 'Optic Radiation FD', 'fdcopticRadiation', 'intracranialVol', 'intracranialHCP', 'TIV', 'LGN', 'meanAdjustedGCVol', 'V1surface', 'V1thickness'};

fprintf('\n<strong>Correlation of each variable\n</strong>')
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
        
        % Make bar plots for correlations
        if strcmp(measureSet{ii}, 'Height_inches')
            ValheightR = R(1,2);
            ValheightConfidence = std(bootR);
        end
        if strcmp(measureSet{ii}, 'Weight_pounds')
            ValweightR = R(1,2);
            ValweightConfidence = std(bootR);  
        end
        if strcmp(measureSet{ii}, 'intracranialVol')
            ValicvR = R(1,2);
            ValicvConfidence = std(bootR);  
        end
    end
    err = [ValheightConfidence ValweightConfidence ValicvConfidence];
    figure;
    X = categorical({'height' 'weight' 'ICV'});
    X = reordercats(X,{'height' 'weight' 'ICV'});
    y = [ValheightR, ValweightR, ValicvR];
    bar(X,y)
    ylim([-0.5 0.5])
    ylabel('R')
    title([names{ff} ' and size correlations'])
    hold on
    er = errorbar(X,y,err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    hold off
    ValheightR = [];
    ValweightR = [];
    ValicvR = [];
end

%% PCA of the biometric variables 
% pcaMat = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.TIV fixelComparisonTable.Axial_Length_average];
% categoryNames = ["height"; "weight"; "ICV"; "AxialLength"];

pcaMat = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.TIV];
categoryNames = ["height"; "weight"; "ICV"];

[wcoeff,score,latent,tsquared,explained] = pca(pcaMat,'VariableWeights','variance');
coefforth = inv(diag(std(pcaMat)))*wcoeff;

I = coefforth'*coefforth;
I(1:3,1:3)

cscores = zscore(pcaMat)*coefforth;

figure()
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

figure()
biplot(coefforth(:,1:2),'Scores',score(:,1:2),'Varlabels',categoryNames);

PC1 = score(:,1);
PC2 = score(:,2);

% %% PCA for brain mesurements
% pcaMat = [fixelComparisonTable.fc_ fixelComparisonTable.LGN fixelComparisonTable.fc_opticRadiation fixelComparisonTable.V1surface fixelComparisonTable.V1thickness];
% categoryNames = ["fc"; "LGN"; "opticRadiationFC"; "V1surface"; "V1thickness"];
% 
% [wcoeff,score,latent,tsquared,explained] = pca(pcaMat,'VariableWeights','variance');
% coefforth = inv(diag(std(pcaMat)))*wcoeff;
% 
% I = coefforth'*coefforth;
% I(1:3,1:3)
% 
% cscores = zscore(pcaMat)*coefforth;
% 
% figure()
% pareto(explained)
% xlabel('Principal Component')
% ylabel('Variance Explained (%)')
% 
% figure()
% biplot(coefforth(:,1:2),'Scores',score(:,1:2),'Varlabels',categoryNames);
% 
% PC1brain = score(:,1);
% PC2brain = score(:,2);
% 
% [r, p] = corrcoef(PC1, PC1brain);
% r = num2str(r(2));
% p = num2str(p(2));
% fprintf(['Correlation of biometric PC1 with the PC1 of the brain measurements r=' r ' ' 'p=' p '\n'])  
% [r, p] = corrcoef(fixelComparisonTable.meanAdjustedGCVol, PC1brain);
% r = num2str(r(2));
% p = num2str(p(2));
% fprintf(['Correlation of meanAdjustedGC with the PC1 of the brain measurements r=' r ' ' 'p=' p])  


%% Report partial correlation results
% % Reshape gender matrix so Males will be one and females will be zero
% genderMatrix = fixelComparisonTable.Gender;
% for ii = 1:length(genderMatrix)
%     if strcmp(genderMatrix{ii}, 'M')
%         genderMatrix{ii} = 1;
%     else
%         genderMatrix{ii} = 0;
%     end
% end
% genderMatrix = cell2mat(genderMatrix);
% 
% % Define size matrix 
% sizeMatrix = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.TIV fixelComparisonTable.Age genderMatrix];

fprintf('\n<strong>Correlation of GC volume with all structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with opticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.V1thickness, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with V1thickness controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Correlation of optic Tract FC with all other structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation optic Tract FC with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation optic Tract FC with opticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation optic Tract FC with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation optic Tract FC with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Correlation of optic Tract FD with all other structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.fd_, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation optic Tract FD with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation optic Tract FD with opticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation optic Tract FD with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation optic Tract FD with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])


fprintf('\n<strong>Correlation of LGN with all other structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.LGN, fixelComparisonTable.fc_, PC1);
fprintf(['\nPartial correlation LGN with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation LGN with opticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.LGN, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation LGN with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.LGN, fixelComparisonTable.V1thickness, PC1);
fprintf(['\nPartial correlation LGN with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Correlation of Optic Radiation FC with all other structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.fc_, PC1);
fprintf(['\nPartial correlation Optic Radiation FC with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation Optic Radiation FC with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation Optic Radiation FC with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation Optic Radiation FC with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Correlation of V1 surface with all other structures</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.V1surface, fixelComparisonTable.fc_, PC1);
fprintf(['\nPartial correlation V1 surface with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.V1surface, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation V1 surface with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.V1surface, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation V1 surface with opticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.V1surface, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation V1 surface with V1surface controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Correlation of adjacent structures FC</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.new_LGN, PC1);
fprintf(['\nPartial correlation opticTractFC with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.new_LGN, fixelComparisonTable.fc_opticRadiation, PC1);
fprintf(['\nPartial correlation LGN with OpticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation OpticRadiationFC with V1surfaceArea controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1thickness, PC1);
fprintf(['\nPartial correlation OpticRadiationFC with V1thickness controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n\n'])
% [rho, pval] = partialcorr(fixelComparisonTable.fc_, fixelComparisonTable.LGN, PC1);
% fprintf(['\nPartial correlation opticTractFC with HCP-LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
% [rho, pval] = partialcorr(fixelComparisonTable.new_LGN, fixelComparisonTable.fc_opticRadiation, PC1);
% fprintf(['\nPartial correlation HCP-LGN with OpticRadiationFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])

fprintf('\n<strong>Correlation of adjacent structures FD</strong>')
[rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fd_, PC1);
fprintf(['\nPartial correlation meanAdjustedGCVol with opticTractFD controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_, fixelComparisonTable.LGN, PC1);
fprintf(['\nPartial correlation opticTractFD with LGN controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.LGN, fixelComparisonTable.fd_opticRadiation, PC1);
fprintf(['\nPartial correlation LGN with OpticRadiationFD controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.V1surface, PC1);
fprintf(['\nPartial correlation OpticRadiationFD with V1surfaceArea controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
[rho, pval] = partialcorr(fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.V1thickness, PC1);
fprintf(['\nPartial correlation OpticRadiationFD with V1thickness controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval) '\n'])

fprintf('\n<strong>Chiasm version</strong>')
% [rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.opticChiasmVolume_Min, PC1);
% fprintf(['\nPartial correlation meanAdjustedGCVol with opticChiasm controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
% [rho, pval] = partialcorr(fixelComparisonTable.opticChiasmVolume_Min, fixelComparisonTable.fc_, PC1);
% fprintf(['\nPartial correlation opticChiasm with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
% [rho, pval] = partialcorr(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.opticChiasmVolume_Min, PC1, 'rows', 'pairwise');
% fprintf(['\nPartial correlation meanAdjustedGCVol with MINopticChiasm controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])
% [rho, pval] = partialcorr(fixelComparisonTable.opticChiasmVolume_Min, fixelComparisonTable.fc_, PC1, 'rows', 'pairwise');
% fprintf(['\nPartial correlation MINopticChiasm with opticTractFC controlled for PC1: ' 'rho:' num2str(rho) ', p:' num2str(pval)])


% %% Remove the effect of biometric PC1 from the brain vectors, then do PCA on the residuals, and then ask if PC1 of the (corrected for bigness) brain vectors is correlated with adjusted GC volume
% 
% fprintf('\n<strong>Brain vectors correlation with adjustedGC after biometric PC1 were partialed out</strong>')
% 
% corrTable = fitlm(fixelComparisonTable.fc_, PC1);
% residualFC = corrTable.Residuals.Pearson;
% 
% corrTable = fitlm(fixelComparisonTable.LGN, PC1);
% residualLGN = corrTable.Residuals.Pearson;
% 
% corrTable = fitlm(fixelComparisonTable.fc_opticRadiation, PC1);
% residualfc_opticRadiation = corrTable.Residuals.Pearson;
% 
% corrTable = fitlm(fixelComparisonTable.V1surface, PC1);
% residualV1surface = corrTable.Residuals.Pearson;
% 
% corrTable = fitlm(fixelComparisonTable.V1thickness, PC1);
% residualV1thickness = corrTable.Residuals.Pearson;
% 
% biometricResidualSet = [residualFC residualLGN residualfc_opticRadiation residualV1surface residualV1thickness];
% categoryNames = ["residualFC"; "residualLGN"; "residualfc_opticRadiation"; "residualV1surface"; "residualV1thickness"];
% 
% [wcoeff,score,latent,tsquared,explained] = pca(biometricResidualSet,'VariableWeights','variance');
% I = coefforth'*coefforth;
% I(1:3,1:3)
% 
% cscores = zscore(pcaMat)*coefforth;
% 
% figure()
% pareto(explained)
% xlabel('Principal Component')
% ylabel('Variance Explained (%)')
% 
% figure()
% biplot(coefforth(:,1:2),'Scores',score(:,1:2),'Varlabels',categoryNames);
% 
% PC1residuals = score(:,1);
% PC2residuals = score(:,2);
% 
% [p, r] = corrcoef(PC1residuals, fixelComparisonTable.meanAdjustedGCVol);
% r = num2str(r(2));
% p = num2str(p(2));
% fprintf(['Correlation of size regressed brainVector PC1 with the adjustedGCVol r=' r ' ' 'p=' p '\n'])  
%% Controlled correlation plots showing the correlation of adjacent regions
x = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.fc_, fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation];
y = [fixelComparisonTable.fc_, fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface];
xName = {'meanAdjustedGCVol', 'OpticTractFC', 'LGN Volume', 'OpticRadiationFC'};
yName = {'OpticTractFC','LGN Volume', 'OpticRadiationFC', 'V1surface'}; 

for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    [R,P] = corrcoef(residualsx, residualsy);
    mycorr = @(residualsx,residualsy) corr(residualsx,residualsy);
    niterations = 10000;
    interval = bootci(niterations,{mycorr,residualsx,residualsy});
    RL = interval(1);
    RU = interval(2);
    
    % add first plot in 2 x 1 grid  
    figure
    scatter(residualsx, residualsy, 'MarkerFaceColor', 'k');
    xlim([-3 3])
    ylim([-3 3])
    xlabel([xName(ii) ' size PC1 partialled out']);
    ylabel([yName(ii) ' size PC1 partialled out']);
    title([xName(ii) ' vs ' yName(ii) ' controlled for height/weight/ICV/Gender/Age'])
    box 'on'
    axis square;
    set(gca,'Ticklength',[0 0])
    %white background
    set(gcf,'color','w');
    refline
    theStringR = sprintf(['R=' ' ' num2str(sprintf('%.3f', R(1,2)))], residualsx,  residualsy);
    theStringP = sprintf(['P=' ' ' num2str(sprintf('%.3f', P(1,2)))], residualsx,  residualsy);
    theStringRL = sprintf(['RL=' ' ' num2str(sprintf('%.3f', RL))], residualsx,  residualsy);
    theStringRU = sprintf(['RU=' ' ' num2str(sprintf('%.3f', RU))], residualsx,  residualsy);
    text(1, -1, theStringR, 'FontSize', 10);
    text(1, -1.4, theStringP, 'FontSize', 10);
    text(1, -1.8, theStringRL, 'FontSize', 10);
    text(1, -2.2, theStringRU, 'FontSize', 10);
end

%% Controlled correlation plots showing the correlation between GC and everything else
x = [fixelComparisonTable.meanAdjustedGCVol,fixelComparisonTable.meanAdjustedGCVol,fixelComparisonTable.meanAdjustedGCVol,fixelComparisonTable.meanAdjustedGCVol];
y = [fixelComparisonTable.fc_, fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface];
xName = {'meanAdjustedGCVol', 'meanAdjustedGCVol', 'meanAdjustedGCVol', 'meanAdjustedGCVol'};
yName = {'OpticTractFC','LGN Volume', 'OpticRadiationFC', 'V1surface'}; 

for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    [R,pval,RL,RU] = corrcoef(residualsx, residualsy);
    figure;
    % add first plot in 2 x 1 grid    
    scatter(residualsx, residualsy, 'MarkerFaceColor', 'k');
    xlim([-3 3])
    ylim([-3 3])
    xlabel(['Relative ' xName(ii) ' to size']);
    ylabel(['Relative ' yName(ii) ' to size']);
    title([xName(ii) ' vs ' yName(ii) ' controlled for height/weight/ICV/Gender/Age'])
    box 'on'
    axis square;
    set(gca,'Ticklength',[0 0])
    %white background
    set(gcf,'color','w');
    refline
    theStringR = sprintf(['R=' ' ' num2str(sprintf('%.3f', R(1,2)))], residualsx,  residualsy);
    theStringP = sprintf(['P=' ' ' num2str(sprintf('%.3f', pval(1,2)))], residualsx,  residualsy);
    theStringRL = sprintf(['RL=' ' ' num2str(sprintf('%.3f', RL(1,2)))], residualsx,  residualsy);
    theStringRU = sprintf(['RU=' ' ' num2str(sprintf('%.3f', RU(1,2)))], residualsx,  residualsy);
    text(1, -1, theStringR, 'FontSize', 10);
    text(1, -1.4, theStringP, 'FontSize', 10);
    text(1, -1.8, theStringRL, 'FontSize', 10);
    text(1, -2.2, theStringRU, 'FontSize', 10);
end

% %% Mediation value 
% x = [fixelComparisonTable.fc_, fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation];
% y = [fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface];
% z = fixelComparisonTable.meanAdjustedGCVol;
% xName = {'OpticTractFC', 'LGNVolume', 'OpticRadiationFC'};
% yName = {'LGNVolume', 'OpticRadiationFC', 'V1surface'}; 
% zName = {'meanAdjustedGCVol','meanAdjustedGCVol','meanAdjustedGCVol'};
% 
% corrTablez = fitlm(sizeMatrix, z);
% residualsz = corrTablez.Residuals.Pearson;
% 
% for ii = 1:length(xName)
%     corrTablex = fitlm(sizeMatrix, x(1:end, ii));
%     residualsx = corrTablex.Residuals.Pearson;
% 
%     corrTabley = fitlm(sizeMatrix, y(1:end, ii));
%     residualsy = corrTabley.Residuals.Pearson;
%     
%     mediation(residualsx, residualsy, residualsz, 'plots', 'names', {xName(ii) yName(ii) zName(ii)}, 'dosave')
%     currentFolder = pwd;
%     pathDiagramPng = fullfile(currentFolder, 'Path_Diagram.png');
%     pathDiagramFig = fullfile(currentFolder, 'Path_Diagram.fig');
%     scatterPng = fullfile(currentFolder, 'Mediation_Scatterplots.png');
%     scatterFig = fullfile(currentFolder, 'Mediation_Scatterplots.fig');
%     textFile = fullfile(currentFolder, 'Mediation_Output.txt');
%     mediationPlots = fullfile(p.Results.fixelDataDir, 'mediationPlots');
%     if ~exist(mediationPlots, 'dir')
%         system(['mkdir' ' ' mediationPlots]);
%     end
%     system(['mv ' pathDiagramPng ' ' fullfile(mediationPlots, [xName{ii} '_vs_' yName{ii} '_PathDiagram' '.png'])]);
%     system(['mv ' pathDiagramFig ' ' fullfile(mediationPlots, [xName{ii} '_vs_' yName{ii} '_PathDiagram' '.fig'])]);
%     system(['mv ' scatterPng ' ' fullfile(mediationPlots, [xName{ii} '_vs_' yName{ii} '_Scatter' '.png'])]);
%     system(['mv ' scatterFig ' ' fullfile(mediationPlots, [xName{ii} '_vs_' yName{ii} '_Scatter' '.fig'])]);
%     system(['mv ' textFile ' ' fullfile(mediationPlots, [xName{ii} '_vs_' yName{ii} '_Output' '.txt'])]);
% end

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

% %% Plot fc by LGN values
% 
% [R,pval] = corrcoef(fixelComparisonTable.fc_, fixelComparisonTable.LGN);
% figure;
% % add first plot in 2 x 1 grid    
% scatter(fixelComparisonTable.fc_, fixelComparisonTable.LGN, 'MarkerFaceColor', 'k');
% xlabel ('FC');
% ylabel('LGN');
% title('FC vs LGN')
% box 'on'
% axis square;
% set(gca,'Ticklength',[0 0])
% %white background
% set(gcf,'color','w');
% refline
% theStringR = sprintf(['R=' ' ' num2str(R(1,2))], fixelComparisonTable.fc_,  fixelComparisonTable.LGN);
% theStringP = sprintf(['P=' ' ' num2str(pval(1,2))], fixelComparisonTable.fc_,  fixelComparisonTable.LGN);
% text(1.5, -2.5, theStringR, 'FontSize', 10);
% text(1.5, -2.8, theStringP, 'FontSize', 10);
% 
% %% Plot fd by LGN values
% 
% [R,pval] = corrcoef(fixelComparisonTable.fd_, fixelComparisonTable.LGN);
% figure;
% % add first plot in 2 x 1 grid    
% scatter(fixelComparisonTable.fd_, fixelComparisonTable.LGN, 'MarkerFaceColor', 'k');
% xlabel ('FD');
% ylabel('LGN');
% title('FD vs LGN')
% box 'on'
% axis square;
% set(gca,'Ticklength',[0 0])
% %white background
% set(gcf,'color','w');
% refline
% theStringR = sprintf(['R=' ' ' num2str(R(1,2))], fixelComparisonTable.fd_,  fixelComparisonTable.LGN);
% theStringP = sprintf(['P=' ' ' num2str(pval(1,2))], fixelComparisonTable.fd_,  fixelComparisonTable.LGN);
% text(1.5, -2.5, theStringR, 'FontSize', 10);
% text(1.5, -2.8, theStringP, 'FontSize', 10);
% 
% %% Plot adjustedGC by LGN values
% 
% [R,pval] = corrcoef(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN);
% figure;
% % add first plot in 2 x 1 grid    
% scatter(fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, 'MarkerFaceColor', 'k');
% xlabel ('AdjustedGCVol');
% ylabel('LGNVol');
% title('AdjustedGCVol and LGN')
% box 'on'
% axis square;
% set(gca,'Ticklength',[0 0])
% %white background
% set(gcf,'color','w');
% refline
% theStringR = sprintf(['R=' ' ' num2str(R(1,2))], fixelComparisonTable.fd_,  fixelComparisonTable.LGN);
% theStringP = sprintf(['P=' ' ' num2str(pval(1,2))], fixelComparisonTable.fd_,  fixelComparisonTable.LGN);
% text(2.1*10^-3, 260, theStringR, 'FontSize', 10);
% text(2.1*10^-3, 250, theStringP, 'FontSize', 10);
% 
% %% Plot FC left right
% 
% [R,pval] = corrcoef(fixelComparisonTable.left_fc_, fixelComparisonTable.right_fc_);
% figure;
% % add first plot in 2 x 1 grid    
% scatter(fixelComparisonTable.left_fc_, fixelComparisonTable.right_fc_, 'MarkerFaceColor', 'k');
% xlabel ('leftHemiFC');
% ylabel('rightHemiFC');
% title('Left vs right hemi optic tract FC')
% box 'on'
% axis square;
% set(gca,'Ticklength',[0 0])
% %white background
% set(gcf,'color','w');
% refline
% theStringR = sprintf(['R=' ' ' num2str(sprintf('%.3f', R(1,2)))], fixelComparisonTable.left_fc_,  fixelComparisonTable.right_fc_);
% text(1.1, 0.82, theStringR, 'FontSize', 10);untitled
% 
% %% Plot FD left right
% 
% [R,pval] = corrcoef(fixelComparisonTable.left_fd_, fixelComparisonTable.right_fd_);
% figure;
% % add first plot in 2 x 1 grid    
% scatter(fixelComparisonTable.left_fd_, fixelComparisonTable.right_fd_, 'MarkerFaceColor', 'k');
% xlabel ('leftHemiFD');
% ylabel('rightHemiFD');
% title('Left vs right hemi optic tract FD')
% box 'on'
% axis square;
% set(gca,'Ticklength',[0 0])
% %white background
% set(gcf,'color','w');
% refline
% theStringR = sprintf(['R=' ' ' num2str(sprintf('%.3f', R(1,2)))], fixelComparisonTable.left_fd_,  fixelComparisonTable.right_fd_);
% text(0.68, 0.52, theStringR, 'FontSize', 10);
% 
% %% Bar plot for correlations 
% barItems = {'fc_', 'fd_', 'meanAdjustedGCVol'};
% barNames = {'FC', 'FD', 'meanAdjustedGCVol'};
% for b = 1:length(barItems)
%     figure;
%     FCheight = corrcoef(fixelComparisonTable.(barItems{b}), fixelComparisonTable.Height_inches);
%     FCheight = FCheight(2);
%     FCweight = corrcoef(fixelComparisonTable.(barItems{b}), fixelComparisonTable.Weight_pounds);
%     FCweight = FCweight(2);
%     FCICV = corrcoef(fixelComparisonTable.(barItems{b}), fixelComparisonTable.intracranialVol); 
%     FCICV = FCICV(2);
%     X = categorical({'height' 'weight' 'ICV'});
%     X = reordercats(X,{'height' 'weight' 'ICV'});
%     y = [FCheight, FCweight, FCICV];
%     bar(X,y)
%     ylim([-0.5 0.5])
%     ylabel('R')
%     title([barNames{b} ' and size correlations'])
% end
end

%% LOCAL FUNCTIONS

function surfaceArea = calcSurfaceArea(vert,face,vertIdx)

% Find the faces that are composed entirely of vertices in the index
faceIdx = find(all(ismember(face,vertIdx)'));

v1 = vert(face(faceIdx,2),:)-vert(face(faceIdx,1),:);
v2 = vert(face(faceIdx,3),:)-vert(face(faceIdx,2),:);
cp = 0.5*cross(v1,v2);
surfaceAreaAll = sum(sqrt(dot(cp,cp,2)));

% Find the faces that are composed of any vertices in the index
faceIdx = find(any(ismember(face,vertIdx)'));

v1 = vert(face(faceIdx,2),:)-vert(face(faceIdx,1),:);
v2 = vert(face(faceIdx,3),:)-vert(face(faceIdx,2),:);
cp = 0.5*cross(v1,v2);
surfaceAreaAny = sum(sqrt(dot(cp,cp,2)));

% Report the average of these two
surfaceArea = (surfaceAreaAll+surfaceAreaAny)/2;

end

