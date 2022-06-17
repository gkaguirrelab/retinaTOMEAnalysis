function fixelAnalysis(varargin)
%% Set the dropboxBaseDir and flywheel id
% We need this for the default loations of some the directories
% dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));
dropboxBaseDir='/home/ozzy/Dropbox (Aguirre-Brainard Lab)';

fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Get tome subjects
projects = fw.projects();
tome = projects{1,1};
subjects = tome.subjects();
subjectLength = length(subjects);

% Set rng seed
rng('default')

% % Download freesurfer zips in case we need them
% freesurfer_subject_path = '/home/ozzy/freesurfer_subjects';
% for sub = 1:subjectLength
%     if strcmp(subjects{sub}.label(1:4), 'TOME') && ~strcmp(subjects{sub}.label(6:end), '3027')      
%         % Get subject name and save the name for where aseg files will be saved
%         subject = subjects{sub,1};
%         subjectLabel = subject.label;
%         % Do the next block if any of the stat files do not exist in the path
%         sessions = subject.sessions();
%         for ses = 1:length(sessions)
%             session = sessions{ses,1};
%             % Get analysis and loop through
%             analyses = session.analyses();
%             for a = 1:length(analyses)
%                 % Get aseg from latest freesurfer runs     
%                 if contains(analyses{a,1}.label, 'freesurfer')
%                     freesurferAnalysisContainer = analyses{a,1};
%                     analysisTag = freesurferAnalysisContainer.id;
%                     zipFile = ['freesurfer-recon-all_' subjectLabel '_' analysisTag '.zip'];  
%                     freesurferAnalysisContainer.downloadFile(zipFile, fullfile(freesurfer_subject_path, [subjectLabel '.zip']))
%                     unzip(fullfile(freesurfer_subject_path, [subjectLabel '.zip']), freesurfer_subject_path)
%                     delete(fullfile(freesurfer_subject_path, [subjectLabel '.zip']))
%                 end
%             end
%         end
%     end
% end                            
    
%% GCvolume corrections
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
XPos_Degs = [XPos_Degs_horiz; XPos_Degs_vert+posShift];

% Conduct PCA upon the tissue volume data
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

%% Download the fixel values for optic tract and radiations from Flywheel
% Right then left optic tract
laterality = {'right','left'};
analysisIDs = {'613222eae4575d3ae7e439e1','613222e90cb137b99763de1d'};
% analysisIDs = {'60a33c76b4a131197e7bfaa8','60a33c7617fcfbb03ffeacf6'};
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
analysisIDs = {'613dd3a1af9c5aae928f7b26','613dd3a150fe1777a6965bad'};
% analysisIDs = {'60ed117b6d2438c15c96c6bd','60ed1153dcf573726496c77d'}; % NOT A SINGLE SHELL
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

%% Optic tract fixel correlations with multi shell
% Right then left optic tract
laterality = {'right','left'};
analysisIDs = {'60a33c76b4a131197e7bfaa8','60a33c7617fcfbb03ffeacf6'};
fileNames = {'fc_stats.csv','fd_stats.csv','fdc_stats.csv'};
for ll = 1:length(laterality)
    for ff = 1:length(fileNames)
        saveName = fullfile(p.Results.fixelDataDir,[laterality{ll} '_multiShell_' fileNames{ff}]);
        fw.downloadOutputFromAnalysis(analysisIDs{ll},fileNames{ff},saveName);
        
        % Now load the file
        opts = detectImportOptions(saveName);
        fixelData = readtable(saveName, opts);
        if ll==1 && ff==1
            fixelTableMultiShell = fixelData(:,1:2);
            fixelTableMultiShell.Properties.VariableNames{2} = [laterality{ll} '_multiShell_' fileNames{ff}(1:3)];
        else
            subTable = fixelData(:,1:2);
            subTable.Properties.VariableNames{2} = [laterality{ll} '_multiShell_' fileNames{ff}(1:3)];
            fixelTableMultiShell=join(fixelTableMultiShell,subTable);
        end
    end
end

% Massage the fixelTable to match up with the comboTable
fixelTableMultiShell.Properties.VariableNames{1} = 'TOME_ID';
fixelTableMultiShell.TOME_ID = strrep(fixelTableMultiShell.TOME_ID,'fod_','');

% Sort rows by subject ID, so that it will be easier to add other measures
fixelTableMultiShell = sortrows(fixelTableMultiShell);
fixelTable = join(fixelTable, fixelTableMultiShell);

%% Download the DTI values for optic tract and radiations from Flywheel
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

laterality = {'right','left','right','left'};
analysisIDs = {'623a32bc3b0ec157fb571815','623a31ab0f25db25cbcc3e17','623a314e76a5a3f9ad74295b','623a3074bc08ba8422ca2906'};
fileNames = {'FA_stats.csv', 'FA_stats.csv', 'MD_stats.csv', 'MD_stats.csv'};
for ll = 1:length(laterality)
    saveName = fullfile(p.Results.fixelDataDir,[laterality{ll} '_OR' fileNames{ll}]);
    fw.downloadOutputFromAnalysis(analysisIDs{ll},fileNames{ll},saveName);
        
    % Now load the file
    opts = detectImportOptions(saveName);
    dtiData = readtable(saveName, opts);
        
    % Sort the table
    dtiData = sortrows(dtiData);
        
    % Add 
    subTable = dtiData(:,1:2);
    subTable.Properties.VariableNames{1} = 'TOME_ID';
    subTable.Properties.VariableNames{2} = [laterality{ll} '_OR' fileNames{ll}(1:2)];
    if strcmp(fixelTable{1,1}{1}(1:7), 'preproc')
        for ii = 1:height(fixelTable)
            fixelTable{ii,1}{1} = fixelTable{ii,1}{1}(12:end);
        end
    end
    fixelTable=join(fixelTable,subTable);
end

%% Get the measurements from all other brain regions and save them into the table
% Create a subject data folder where the output will be saved
subjectDataFolder = fullfile(p.Results.fixelDataDir, 'subjectDataFolder');
if ~exist(subjectDataFolder, 'dir')
    system(['mkdir' ' ' subjectDataFolder]);
end

% Create empty matrices for each of our measurements
subjectNames = [];
leftLGN = [];
rightLGN = [];
V1surfaceLeft = [];
V1surfaceRight = [];
V1ThicknessLeft = [];
V1ThicknessRight = [];
radiationTable = [];
totalSurfaceArea = [];

% Loop through subjects
for sub = 1:subjectLength
    % If label starts with TOME and it's not TOME_3027 (because that
    % subject was discarted), process the subject
    if strcmp(subjects{sub}.label(1:4), 'TOME') && ~strcmp(subjects{sub}.label(6:end), '3027')      
        % Get subject name and save the name for where aseg files will be saved
        subject = subjects{sub,1};
        subjectLabel = subject.label;
        subjectFolder = fullfile(subjectDataFolder, subjectLabel);
        if ~exist(subjectFolder, 'dir')
            system(['mkdir' ' ' subjectFolder]);
        end
        subjectNames = [subjectNames; {subjectLabel}];
        asegSaveName = fullfile(subjectFolder,[subjectLabel '_FSaseg.stats']);
        asegSaveNameHCP = fullfile(subjectFolder,[subjectLabel '_HCPaseg.stats']);
        lhWhite = fullfile(subjectFolder,[subjectLabel '_lh.white']);
        rhWhite = fullfile(subjectFolder,[subjectLabel '_rh.white']);        
        lhThickness = fullfile(subjectFolder,[subjectLabel '_lh.thickness']);
        rhThickness = fullfile(subjectFolder,[subjectLabel '_rh.thickness']);
        LGNSaveName = fullfile(subjectFolder,[subjectLabel '_FSThalamicNuclei.v12.T1.volumes.txt']);
        bayesPrfSaveNameLeft = fullfile(subjectFolder,['lh.' subjectLabel '_inferred_varea.mgz']);
        bayesPrfSaveNameRight = fullfile(subjectFolder,['rh.' subjectLabel '_inferred_varea.mgz']);
        opticRadiationMaskSaveName = fullfile(subjectFolder,[subjectLabel '_mask_combined_in_FOD_template.nii.gz']);

        % Do the next block if any of the stat files do not exist in the path
        if ~isfile(asegSaveName) || ~isfile(asegSaveNameHCP) || ~isfile(lhWhite) || ~isfile(rhWhite) ||~isfile(rhWhite) || ~isfile(lhThickness) || ~isfile(rhThickness) || ~isfile(LGNSaveName) || ~isfile(bayesPrfSaveNameLeft) || ~isfile(bayesPrfSaveNameRight) || ~isfile(opticRadiationMaskSaveName)
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
                    % get thickness, white surfaces, and aseg
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
                    % Find analyses that contain bayesPRF in the name
                    if ~isfile(bayesPrfSaveNameLeft)
                        if contains(analyses{a,1}.label, 'bayesprf')
                            freesurferAnalysisContainer = analyses{a,1};
                            zipFile = [subjectLabel '_inferred_surface.zip'];
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, ['lh.' subjectLabel '_inferred_varea.mgz'], bayesPrfSaveNameLeft);
                            freesurferAnalysisContainer.downloadFileZipMember(zipFile, ['rh.' subjectLabel '_inferred_varea.mgz'], bayesPrfSaveNameRight);
                        end
                    end
                    % Download the opticRadiation volumetric masks 
                    if ~isfile(opticRadiationMaskSaveName)
                        if contains(analyses{a,1}.label, 'SS3T - traceOpticRadiation')
                            freesurferAnalysisContainer = analyses{a,1};
                            fileName = [subjectLabel '_mask_combined_in_FOD_template.nii.gz'];
                            file = freesurferAnalysisContainer.getFile(fileName);
                            file.download(opticRadiationMaskSaveName);
                        end
                    end
                end
            end
        end               
        
        % Load FS Thalamic segmentation text and extract LGN values
        fid = fopen(LGNSaveName);
        leftLGNLine = fgetl(fid);
        leftLGNVal = str2num(leftLGNLine(10:end));
        rightLGNLine = fgetl(fid);
        rightLGNVal = str2num(rightLGNLine(11:end));
        leftLGN = [leftLGN; leftLGNVal];
        rightLGN = [rightLGN; rightLGNVal];
        fclose(fid);
            
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
    
        % Total surface area
        [~, outlh] = system(['mris_anatomical_stats ' subjectLabel ' lh']);
        [~, outrh] = system(['mris_anatomical_stats ' subjectLabel ' rh']);
        outlh = extractBetween(outlh, 'total surface area                      = ', ' mm^2');
        outrh = extractBetween(outrh, 'total surface area                      = ', ' mm^2');
        outlh = str2num(outlh{1});
        outrh = str2num(outrh{1});
        totalhemisarea = outlh + outrh;
        totalSurfaceArea = [totalSurfaceArea; totalhemisarea];
       
        % Construct optic radiation volume from 
        radiationMask = MRIread(opticRadiationMaskSaveName);
        lengthNonZero = find(radiationMask.vol);
        lengthNonZero = length(lengthNonZero);
        radiationVolume = lengthNonZero * 1.25;   
        radiationTable = [radiationTable; radiationVolume];
    end
end

% Combine subject and LGN in a table and sort rows.
LGN = (leftLGN + rightLGN) / 2;
LGNTable = table(subjectNames, leftLGN, rightLGN, LGN);
LGNTable = sortrows(LGNTable);
LGNTable.Properties.VariableNames{1} = 'TOME_ID';
LGNTable.Properties.VariableNames{2} = 'left_LGN';
LGNTable.Properties.VariableNames{3} = 'right_LGN';
fixelTable=join(fixelTable,LGNTable);

% % Add total surface to the fixel table
% totalSurfaceArea = table(subjectNames, totalSurfaceArea);
% totalSurfaceArea.Properties.VariableNames{1} = 'TOME_ID';
% totalSurfaceArea.Properties.VariableNames{2} = 'totalSurfaceArea';
% fixelTable=join(fixelTable,totalSurfaceArea);

% Add optic radiation volume to the table
radiationVolumeTable = table(subjectNames, radiationTable);
radiationVolumeTable = sortrows(radiationVolumeTable);
radiationVolumeTable.Properties.VariableNames{1} = 'TOME_ID';
radiationVolumeTable.Properties.VariableNames{2} = 'opticRadiationVolume';
fixelTable=join(fixelTable,radiationVolumeTable);

% Combine subject and V1 volume and thickness in a table
V1surface = (V1surfaceLeft + V1surfaceRight) / 2;
V1SurfaceTable = table(subjectNames, V1surfaceLeft,  V1surfaceRight, V1surface);
V1SurfaceTable = sortrows(V1SurfaceTable);
V1SurfaceTable.Properties.VariableNames{1} = 'TOME_ID';
V1SurfaceTable.Properties.VariableNames{2} = 'left_V1surface';
V1SurfaceTable.Properties.VariableNames{3} = 'right_V1surface';
fixelTable=join(fixelTable,V1SurfaceTable);
V1thickness = (V1ThicknessLeft + V1ThicknessRight) / 2;
V1ThicknessTable = table(subjectNames, V1ThicknessLeft,  V1ThicknessRight, V1thickness);
V1ThicknessTable = sortrows(V1ThicknessTable);
V1ThicknessTable.Properties.VariableNames{1} = 'TOME_ID';
V1ThicknessTable.Properties.VariableNames{2} = 'left_V1thickness';
V1ThicknessTable.Properties.VariableNames{3} = 'right_V1thickness';
fixelTable=join(fixelTable,V1ThicknessTable);

% Get SPM intracranial volume
SPMTIV = fullfile(p.Results.fixelDataDir, 'TIV.txt');
load(SPMTIV)
subSorted = sortrows(subjectNames);
TIVtable = table(subSorted, TIV);
TIVtable.Properties.VariableNames{1} = 'TOME_ID';
fixelTable=join(fixelTable,TIVtable);

%% Create variables for left-right mean and create fixel comparison table
fixelSet = {'fc_','fd_','fdc', 'multiShell_fc_','multiShell_fd_','multiShell_fdc', 'FA', 'MD', 'ORFA', 'ORMD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation', 'LGN', 'V1surface', 'V1thickness'};
for ff = 1:length(fixelSet)
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    fixelTable.(fixelSet{ff}) = mean([fixelValR, fixelValL],2);
end
 
fixelComparisonTable = join(comboTable(ismember(comboTable.TOME_ID,fixelTable.TOME_ID),:),fixelTable,'Keys','TOME_ID');

%% Run PCA analysis on biometric variables 
% Print
fprintf('\n<strong>Correlation of body size with RGC and optic tract fixel measurements\n</strong>')

% Set matrix
pcaMat = [fixelComparisonTable.Height_inches fixelComparisonTable.Weight_pounds fixelComparisonTable.TIV fixelComparisonTable.Age];
categoryNames = ["height"; "weight"; "ICV"; "Age"];

% Run PCA
[wcoeff,score,latent,tsquared,explained] = pca(pcaMat,'VariableWeights','variance', 'Centered', true);
coefforth = inv(diag(std(pcaMat)))*wcoeff;
I = coefforth'*coefforth;
I(1:3,1:3)
cscores = zscore(pcaMat)*coefforth;

% Create a parieto figure
figure()
pareto(explained)
xlabel('Components')
ylabel('Variance Explained (%)')

% Create a biplot 
figure()
biplot(coefforth(:,1:2),'Scores',score(:,1:2),'Varlabels',categoryNames,'Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1]);

% Save the PC1 and PC2 components
PC1 = score(:,1);
PC2 = score(:,2);

% Plot the first component with the fit
mu = mean(pcaMat);
sd = std(pcaMat);
Xz = bsxfun(@minus, pcaMat ,mu)./sd; 
fig = figure;
fig.Renderer='Painters';
subplot(2,2,1);
scatter3(Xz(:,1),Xz(:,2),Xz(:,3),'bo','Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1])
hold on 
quiver3(0,0,0, coefforth(1,1),coefforth(2,1),coefforth(3,1),3, 'k', 'ShowArrowHead',false);
quiver3(0,0,0, -coefforth(1,1),-coefforth(2,1),-coefforth(3,1),3, 'k', 'ShowArrowHead',false);
xlabel('Height [inches]'); ylabel('Weight [lbs]'); zlabel('ICV [cm^3]');
box off 
grid on
hold off

% Correlation of optic tract FC with the PC1
FCopt = fixelComparisonTable.fc_;
fcFit = fitlm(PC1, FCopt);
R = corr(PC1,FCopt);
mycorr = @(PC1,FCopt) corr(PC1,FCopt);
interval = bootci(10000,{mycorr,PC1,FCopt});
fprintf(['correlation Optic tract FC with PC1 is: ' num2str(R) ' ' '95per CI= ' num2str(interval(1)) ' ' num2str(interval(2)) '\n']);

subplot(2,2,3); 
plot(fcFit, 'Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1])
xlabel(['Body Size']);
ylabel(['Optic Tract FC'])
box off
grid on
title('')
theStringR = sprintf(['N=42, R=' num2str(sprintf('%.2f', R))], PC1,  FCopt);
text(-2.8, 1.135, theStringR, 'FontSize', 10);

% Correlation of optic tract FD with the PC1
FDopt = fixelComparisonTable.fd_;
fdFit = fitlm(PC1, FDopt);
R = corr(PC1,FDopt);
mycorr = @(PC1,FDopt) corr(PC1,FDopt);
interval = bootci(10000,{mycorr,PC1,FDopt});
fprintf(['correlation Optic tract FD with PC1 is: ' num2str(R) ' ' '95per CI= ' num2str(interval(1)) ' ' num2str(interval(2)) '\n']);

subplot(2,2,4);
plot(fdFit, 'Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1])
xlabel(['Body Size']);
ylabel(['Optic Tract FD'])
title('')
box off
grid on
theStringR = sprintf(['N=42, R=' num2str(sprintf('%.2f', R))], PC1,  FDopt);
text(-2.8, 0.69, theStringR, 'FontSize', 10);

h=gcf;
set(h,'PaperPositionMode','auto');         
set(h,'PaperOrientation','landscape');
set(h,'Position',[50 50 1200 800]);

% Correlation of RGC with PC1
GCopt = fixelComparisonTable.meanAdjustedGCVol;
R = corr(PC1,GCopt);
mycorr = @(PC1,GCopt) corr(PC1,GCopt);
interval = bootci(10000,{mycorr,PC1,GCopt});
fprintf(['correlation adjusted GC with PC1 is: ' num2str(R) ' ' '95per CI= ' num2str(interval(1)) ' ' num2str(interval(2)) '\n']);

%% Correlation of left right hemisphere values
fprintf('\n<strong>Correlation of left right hemisphere measurements of diffusion variables\n</strong>')
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'ORFA', 'ORMD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation', 'LGN', 'V1surface', 'V1thickness'};
for ff = 1:length(fixelSet)
    % Report the correlation of left and right
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    X = fitlm(PC1, fixelValR);
    Y = fitlm(PC1, fixelValL);
    resx = X.Residuals.Pearson;
    resy = Y.Residuals.Pearson;
    model = fitlm(resx, resy);
    modelcorr = @(resx,resy) corr(resx,resy);
    niterations = 10000;
    interval = bootci(niterations,{modelcorr,resx,resy});
    R = model.Coefficients{2,1};
    fprintf([fixelSet{ff} ' correlation left with right: ' num2str(R) ' ' '95per CI= ' num2str(interval(1)) ' ' num2str(interval(2)) '\n']);
    % Assign the correlations here  to variables after printing we'll need
    % those for calculating correlation ceilings later
    if strcmp(fixelSet{ff}, 'fc_')
        fcLeftRight = R;
    elseif strcmp(fixelSet{ff}, 'fd_')
        fdLeftRight = R;
    elseif strcmp(fixelSet{ff}, 'fdc')
        fdcLeftRight = R;   
    elseif strcmp(fixelSet{ff}, 'fc_opticRadiation')
        fcRadiationLeftRight = R;
    elseif strcmp(fixelSet{ff}, 'fd_opticRadiation')
        fdRadiationLeftRight = R;  
    elseif strcmp(fixelSet{ff}, 'fdcopticRadiation')
        fdcRadiationLeftRight = R;         
    elseif strcmp(fixelSet{ff}, 'LGN')
        LGNLeftRight = R; 
    elseif strcmp(fixelSet{ff}, 'V1surface')
        V1LeftRight = R;
    elseif strcmp(fixelSet{ff}, 'FA')
        FALeftRight = R;
    elseif strcmp(fixelSet{ff}, 'MD')
        MDLeftRight = R;
    elseif strcmp(fixelSet{ff}, 'ORFA')
        ORFALeftRight = R;
    elseif strcmp(fixelSet{ff}, 'ORMD')
        ORMDLeftRight = R;
    end
end

% We know GC left right eye correlation from Min study 
GCLeftRight = 0.64;

%% Controlled correlation plots showing all correlations
fprintf('\n<strong>Correlation of all structures and fixel values\n</strong>')
x = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol, ...
     fixelComparisonTable.fc_, fixelComparisonTable.fc_, fixelComparisonTable.fc_, ...
     fixelComparisonTable.fd_, fixelComparisonTable.fd_, fixelComparisonTable.fd_, ...
     fixelComparisonTable.fdc, fixelComparisonTable.fdc, fixelComparisonTable.fdc, ...
     fixelComparisonTable.LGN, fixelComparisonTable.LGN, fixelComparisonTable.LGN, fixelComparisonTable.LGN, ...
     fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.fdcopticRadiation, ...
     fixelComparisonTable.meanAdjustedGCVol, ...
     fixelComparisonTable.fc_, fixelComparisonTable.fc_opticRadiation];
y = [fixelComparisonTable.fc_, fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, fixelComparisonTable.fd_, fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.fdc, fixelComparisonTable.fdcopticRadiation, ...
     fixelComparisonTable.LGN, fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.LGN, fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.LGN, fixelComparisonTable.fdcopticRadiation, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.fc_opticRadiation, fixelComparisonTable.V1surface, fixelComparisonTable.fd_opticRadiation, fixelComparisonTable.fdcopticRadiation, ...
     fixelComparisonTable.V1surface, fixelComparisonTable.V1surface, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.opticRadiationVolume, ...
     fixelComparisonTable.fd_, fixelComparisonTable.fd_opticRadiation];
 
xName = {'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', ...
         'Optic Tract FC', 'Optic Tract FC', 'Optic Tract FC', ...
         'Optic Tract FD', 'Optic Tract FD', 'Optic Tract FD', ...
         'Optic Tract FDC', 'Optic Tract FDC', 'Optic Tract FDC', ...
         'LGN Volume', 'LGN Volume', 'LGN Volume', 'LGN Volume', ...
         'Optic Radiation FC', 'Optic Radiation FD', 'Optic Radiation FDC', ...
         'Mean Adjusted GC Volume', ...
         'Optic Tract FC', 'Optic Radiation FC'};
yName = {'Optic Tract FC', 'LGN Volume', 'Optic Radiation FC', 'V1surface', 'Optic Tract FD', 'Optic Radiation FD', 'Optic Tract FDC', 'Optic Radiation FDC', ...
         'LGN Volume', 'Optic Radiation FC', 'V1surface', ...
         'LGN Volume', 'Optic Radiation FD', 'V1surface', ...
         'LGN Volume', 'Optic Radiation FDC', 'V1surface', ...
         'Optic Radiation FC', 'V1surface', 'Optic Radiation FD', 'Optic Radiation FDC', ...
         'V1surface', 'V1surface', 'V1surface', ...
         'Optic Radiation Volume', ...
         'Optic Tract FD', 'Optic Radiation FD'};

fig2 = figure;
fig2.Renderer='Painters';    
subplotNum = 1;
for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    mdl = fitlm(residualsx, residualsy);
    mycorr = @(residualsx,residualsy) corr(residualsx,residualsy);
    niterations = 10000;
    interval = bootci(niterations,{mycorr,residualsx,residualsy});
    R = corr(residualsx, residualsy);
    RL = interval(1);
    RU = interval(2);
    fprintf(['correlation of ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);

    if subplotNum < 5
        % add first plot in 2 x 1 grid  
        subplot(2,2,subplotNum)
        plot(mdl, 'Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1])
        xlabel([xName(ii)]);
        ylabel([yName(ii)]);  
        title('')
        legend off
        box off
        grid on
        theStringR = sprintf(['N=42, R=' num2str(sprintf('%.2f', mdl.Coefficients{2,1}))], residualsx,  residualsy);
        text(-2.8, 2.5, theStringR, 'FontSize', 10);
        xlim([-3 3])
        ylim([-3 3])

        subplotNum = subplotNum + 1;
    end
end

h=gcf;
set(h,'PaperPositionMode','auto');

%% Plot the correlations for not axial corrected GC, but still body size corrected
fprintf('\n<strong>Correelations for non corrected RGC\n</strong>')

% GC not axial corrected, but GC,FC all body corrected
bodyCorrectedFC = fitlm(PC1, fixelComparisonTable.fc_);
bodyCorrectedFC = bodyCorrectedFC.Residuals.Pearson;
bodyCorrectedFD = fitlm(PC1, fixelComparisonTable.fd_);
bodyCorrectedFD = bodyCorrectedFD.Residuals.Pearson;
bodyCorrectedNonAxialGC = fitlm(PC1, fixelComparisonTable.meanFitGCVol);
bodyCorrectedNonAxialGC = bodyCorrectedNonAxialGC.Residuals.Pearson;

mycorr = @(bodyCorrectedNonAxialGC,bodyCorrectedFC) corr(bodyCorrectedNonAxialGC,bodyCorrectedFC);
niterations = 10000;
interval = bootci(niterations,{mycorr,bodyCorrectedNonAxialGC,bodyCorrectedFC});
R = corr(bodyCorrectedNonAxialGC,bodyCorrectedFC);
RL = interval(1);
RU = interval(2);
fprintf(['correlation of uncorrected RGC and FC (both body size corrected) is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);

mycorr = @(bodyCorrectedNonAxialGC,bodyCorrectedFD) corr(bodyCorrectedNonAxialGC,bodyCorrectedFD);
niterations = 10000;
interval = bootci(niterations,{mycorr,bodyCorrectedNonAxialGC,bodyCorrectedFD});
R = corr(bodyCorrectedNonAxialGC,bodyCorrectedFD);
RL = interval(1);
RU = interval(2);
fprintf(['correlation of uncorrected RGC and FD (both body size corrected) is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);

%% DTI correlations with body size correction
fprintf('\n<strong>DTI correlations after body size adjustment\n</strong>')
x = [fixelComparisonTable.FA, fixelComparisonTable.FA, fixelComparisonTable.FA, fixelComparisonTable.FA, ...
     fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, ...
     fixelComparisonTable.MD, fixelComparisonTable.MD, fixelComparisonTable.MD, fixelComparisonTable.MD, ...
     fixelComparisonTable.ORMD, fixelComparisonTable.ORMD, fixelComparisonTable.ORMD, fixelComparisonTable.ORMD];
y = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, fixelComparisonTable.ORFA, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.FA, fixelComparisonTable.LGN, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, fixelComparisonTable.ORMD, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.MD, fixelComparisonTable.LGN, fixelComparisonTable.V1surface];

xName = {'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', ...
         'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', ...
         'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', ...
         'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD'};

yName = {'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation FA', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract FA', 'LGN Volume', 'V1surface', ...
         'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation MD', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract MD', 'LGN Volume', 'V1surface'};
     
fig2 = figure;
fig2.Renderer='Painters';    
subplotNum = 1;
for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    mdl = fitlm(residualsx, residualsy);
    mycorr = @(residualsx,residualsy) corr(residualsx,residualsy);
    niterations = 10000;
    interval = bootci(niterations,{mycorr,residualsx,residualsy});
    R = corr(residualsx, residualsy);
    RL = interval(1);
    RU = interval(2);
    fprintf(['correlation of ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);

    if subplotNum < 5
        % add first plot in 2 x 1 grid  
        subplot(2,2,subplotNum)
        plot(mdl, 'Marker', 'o', 'MarkerEdgeColor',[0, 0, 1], 'MarkerFaceColor',[0, 0, 1])
        xlabel([xName(ii)]);
        ylabel([yName(ii)]);  
        title('')
        legend off
        box off
        grid on
        theStringR = sprintf(['N=42, R=' num2str(sprintf('%.2f', mdl.Coefficients{2,1}))], residualsx,  residualsy);
        text(-2.8, 2.5, theStringR, 'FontSize', 10);
        xlim([-3 3])
        ylim([-3 3])

        subplotNum = subplotNum + 1;
    end
end

h=gcf;
set(h,'PaperPositionMode','auto');     

%% Report DTI correlations without body size correction
fprintf('\n<strong>DTI correlations before the body size correction\n</strong>')
x = [fixelComparisonTable.FA, fixelComparisonTable.FA, fixelComparisonTable.FA, fixelComparisonTable.FA, ...
     fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, fixelComparisonTable.ORFA, ...
     fixelComparisonTable.MD, fixelComparisonTable.MD, fixelComparisonTable.MD, fixelComparisonTable.MD, ...
     fixelComparisonTable.ORMD, fixelComparisonTable.ORMD, fixelComparisonTable.ORMD, fixelComparisonTable.ORMD];
y = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, fixelComparisonTable.ORFA, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.FA, fixelComparisonTable.LGN, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.LGN, fixelComparisonTable.ORMD, fixelComparisonTable.V1surface, ...
     fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.MD, fixelComparisonTable.LGN, fixelComparisonTable.V1surface];

xName = {'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', ...
         'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', ...
         'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', ...
         'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD'};

yName = {'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation FA', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract FA', 'LGN Volume', 'V1surface', ...
         'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation MD', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract MD', 'LGN Volume', 'V1surface'};
     
for ii = 1:length(xName)
    niterations = 10000;
    xCom = x(1:end,ii);
    yCom = y(1:end,ii);
    mycorr = @(xCom,yCom) corr(xCom,yCom);
    R = corr(xCom, yCom);
    interval = bootci(niterations,{mycorr,xCom,yCom});
    RL = interval(1);
    RU = interval(2);
    fprintf(['correlation of ' xName{ii} ' and ' yName{ii}  ' without body size correction is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);
end

h=gcf;
set(h,'PaperPositionMode','auto');

%% Correlation ceilings for fixel data
fprintf('\n<strong>Correlation ceiling values for fixel data\n</strong>')
x = [GCLeftRight, GCLeftRight, GCLeftRight, GCLeftRight, GCLeftRight, GCLeftRight, GCLeftRight, GCLeftRight, ...
     fcLeftRight, fcLeftRight, fcLeftRight, ...
     fdLeftRight, fdLeftRight, fdLeftRight, ...
     fdcLeftRight, fdcLeftRight, fdcLeftRight, ...
     LGNLeftRight, LGNLeftRight, LGNLeftRight, LGNLeftRight, ...
     fcRadiationLeftRight, fdRadiationLeftRight, fdcRadiationLeftRight];
y = [fcLeftRight, LGNLeftRight, fcRadiationLeftRight, V1LeftRight, fdLeftRight, fdRadiationLeftRight, fdcLeftRight, fdcRadiationLeftRight, ...
     LGNLeftRight, fcRadiationLeftRight, V1LeftRight, ...
     LGNLeftRight, fdRadiationLeftRight, V1LeftRight, ...
     LGNLeftRight, fdcRadiationLeftRight, V1LeftRight, ...
     fcRadiationLeftRight, V1LeftRight, fdRadiationLeftRight, fdcRadiationLeftRight, ...
     V1LeftRight, V1LeftRight, V1LeftRight];

xName = {'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume', ...
         'Optic Tract FC', 'Optic Tract FC', 'Optic Tract FC', ...
         'Optic Tract FD', 'Optic Tract FD', 'Optic Tract FD', ...
         'Optic Tract FDC', 'Optic Tract FDC', 'Optic Tract FDC', ...
         'LGN Volume', 'LGN Volume', 'LGN Volume', 'LGN Volume', ...
         'Optic Radiation FC', 'Optic Radiation FD', 'Optic Radiation FDC'};
yName = {'Optic Tract FC', 'LGN Volume', 'Optic Radiation FC', 'V1surface', 'Optic Tract FD', 'Optic Radiation FD', 'Optic Tract FDC', 'Optic Radiation FDC', ...
         'LGN Volume', 'Optic Radiation FC', 'V1surface', ...
         'LGN Volume', 'Optic Radiation FD', 'V1surface', ...
         'LGN Volume', 'Optic Radiation FDC', 'V1surface', ...
         'Optic Radiation FC', 'V1surface', 'Optic Radiation FD', 'Optic Radiation FDC', ...
         'V1surface', 'V1surface', 'V1surface'};
 
for ii = 1:length(x)
    reliabilityA = (x(ii) * 2) / ( x(ii) + 1);
    reliabilityB = (y(ii) * 2) / ( y(ii) + 1);

    corrCeiling = sqrt(reliabilityA*reliabilityB);
    fprintf(['Ceiling values for correlations ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(corrCeiling) '\n']);  
end

%% Correlation ceilings for DTI
fprintf('\n<strong>Correlation ceiling values for dti data\n</strong>')
x = [FALeftRight, FALeftRight, FALeftRight, FALeftRight, ...
     ORFALeftRight, ORFALeftRight, ORFALeftRight, ORFALeftRight, ...
     MDLeftRight, MDLeftRight, MDLeftRight, MDLeftRight, ...
     ORMDLeftRight, ORMDLeftRight, ORMDLeftRight, ORMDLeftRight];
y = [GCLeftRight, LGNLeftRight, ORFALeftRight, V1LeftRight, ...
     GCLeftRight, FALeftRight, LGNLeftRight, V1LeftRight, ...
     GCLeftRight, LGNLeftRight, ORMDLeftRight, V1LeftRight, ...
     GCLeftRight, MDLeftRight, LGNLeftRight, V1LeftRight];

xName = {'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', 'Optic Tract FA', ...
         'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', 'Optic Radiation FA', ...
         'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', 'Optic Tract MD', ...
         'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD', 'Optic Radiation MD'}; 
yName = {'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation FA', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract FA', 'LGN Volume', 'V1surface', ...
         'Mean Adjusted GC Volume', 'LGN Volume', 'Optic Radiation MD', 'V1surface', ...
         'Mean Adjusted GC Volume', 'Optic Tract MD', 'LGN Volume', 'V1surface'};      
     
for ii = 1:length(x)
    reliabilityA = (x(ii) * 2) / ( x(ii) + 1);
    reliabilityB = (y(ii) * 2) / ( y(ii) + 1);

    corrCeiling = sqrt(reliabilityA*reliabilityB);
    fprintf(['Ceiling values for correlations ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(corrCeiling) '\n']);
    
end   

%% Multishell correlations
fprintf('\n<strong>Multi shell correlations\n</strong>')
x = [fixelComparisonTable.meanAdjustedGCVol, fixelComparisonTable.meanAdjustedGCVol];
y = [fixelComparisonTable.multiShell_fc_, fixelComparisonTable.multiShell_fd_];
 
xName = {'Mean Adjusted GC Volume', 'Mean Adjusted GC Volume'};
yName = {'Multi shell FC', 'Multi Shell FD'};

for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    mdl = fitlm(residualsx, residualsy);
    mycorr = @(residualsx,residualsy) corr(residualsx,residualsy);
    niterations = 10000;
    interval = bootci(niterations,{mycorr,residualsx,residualsy});
    R = corr(residualsx, residualsy);
    RL = interval(1);
    RU = interval(2);
    fprintf(['correlation of ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);
end

%% V1 as a  proportion of surface area
V1propSurf = fixelComparisonTable.V1surface./fixelComparisonTable.totalSurfaceArea;
x = [fixelComparisonTable.meanAdjustedGCVol];
y = [V1propSurf];

xName = {'Mean Adjusted GC Volume'};
yName = {'V1/TotalSurface'};

for ii = 1:length(xName)
    corrTablex = fitlm(PC1, x(1:end, ii));
    residualsx = corrTablex.Residuals.Pearson;

    corrTabley = fitlm(PC1, y(1:end, ii));
    residualsy = corrTabley.Residuals.Pearson;

    mdl = fitlm(residualsx, residualsy);
    mycorr = @(residualsx,residualsy) corr(residualsx,residualsy);
    niterations = 10000;
    interval = bootci(niterations,{mycorr,residualsx,residualsy});
    R = corr(residualsx, residualsy);
    RL = interval(1);
    RU = interval(2);
    fprintf(['correlation of ' xName{ii} ' and ' yName{ii}  ' is: ' num2str(R) ' ' '95per CI= ' num2str(RL) ' ' num2str(RU) '\n']);
end
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