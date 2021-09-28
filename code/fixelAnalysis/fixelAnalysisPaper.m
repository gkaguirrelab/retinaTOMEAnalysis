function fixelAnalysisPaper(varargin)

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
XPos_Degs = [XPos_Degs_horiz; XPos_Degs_vert+posShift];

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

%% Get the fixel values for optic tract and chiasm 
% Right then left optic tract
laterality = {'right','left'};
analysisIDs = {'613222eae4575d3ae7e439e1','613222e90cb137b99763de1d'};
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

%% Get the DTI values 
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

%% Get the volume measurements from visual regions and load them to the fixel comparison table
% Create a folder in fixelData dir for extracting aseg files from
% Freesurfer directories
subjectDataFolder = fullfile(p.Results.fixelDataDir, 'subjectDataFolder');
if ~exist(subjectDataFolder, 'dir')
    system(['mkdir' ' ' subjectDataFolder]);
end

% Create empty matrices for everything
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
        new_LGNSaveName = fullfile(subjectFolder,[subjectLabel '_HCPThalamicNuclei.v12.T1.volumes.txt']);
        bayesPrfSaveNameLeft = fullfile(subjectFolder,['lh.' subjectLabel '_inferred_varea.mgz']);
        bayesPrfSaveNameRight = fullfile(subjectFolder,['rh.' subjectLabel '_inferred_varea.mgz']);

        % Do the next block if any of the stat files do not exist in the path
        if ~isfile(asegSaveName) || ~isfile(asegSaveNameHCP) || ~isfile(lhWhite) || ~isfile(rhWhite) ||~isfile(rhWhite) || ~isfile(lhThickness) || ~isfile(rhThickness) || ~isfile(LGNSaveName) || ~isfile(new_LGNSaveName) || ~isfile(bayesPrfSaveNameLeft) || ~isfile(bayesPrfSaveNameRight)
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
        
        % Load FS Thalamic segmentation text and extract LGN values
        fid = fopen(LGNSaveName);
        leftLGNLine = fgetl(fid);
        leftLGNVal = str2num(leftLGNLine(10:end));
        rightLGNLine = fgetl(fid);
        rightLGNVal = str2num(rightLGNLine(11:end));
        leftLGN = [leftLGN; leftLGNVal];
        rightLGN = [rightLGN; rightLGNVal];
        fclose(fid);
       
        % Load LGN Thalamic segmentation text and extract LGN values
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
V1SurfaceTable = table(subjectNames, V1surfaceLeft,  V1surfaceRight, V1surface);
V1SurfaceTable = sortrows(V1SurfaceTable);
V1SurfaceTable.Properties.VariableNames{1} = 'TOME_ID';
V1SurfaceTable.Properties.VariableNames{2} = 'left_V1surface';
V1SurfaceTable.Properties.VariableNames{3} = 'right_V1surface';
fixelTable=join(fixelTable,V1SurfaceTable);

% Combine subject and V1 thickness in a table and sort rows.
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

%% Correlation of left right values
fprintf('\n<strong>Correlation of left right hemisphere measurements of diffusion variables\n</strong>')
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation', 'LGN_HCP', 'V1surface', 'V1thickness'};
for ff = 1:length(fixelSet)
    % Report the correlation of left and right
    fixelValR = fixelTable.(['right_' fixelSet{ff}]);
    fixelValL = fixelTable.(['left_' fixelSet{ff}]);
    R = corr(fixelValR,fixelValL);
    hemicorr = @(fixelValR, fixelValL) corr(fixelValR,fixelValL);
    iterations = 50000;
    interval = bootci(iterations,{hemicorr,fixelValR,fixelValL});  
    fprintf([fixelSet{ff} ' correlation left with right: ' num2str(R) ' ' '95per CI= ' num2str(interval(1)) ' ' num2str(interval(2)) '\n']);
    fixelTable.(fixelSet{ff}) = mean([fixelValR, fixelValL],2);
end

%% Create fixel comparison table 
fixelComparisonTable = join(comboTable(ismember(comboTable.TOME_ID,fixelTable.TOME_ID),:),fixelTable,'Keys','TOME_ID');

%% Make bar plots for biametric values
fixelSet = {'fc_','fd_','fdc', 'FA', 'MD', 'fc_opticRadiation', 'fd_opticRadiation', 'fdcopticRadiation', 'TIV', 'LGN', 'meanAdjustedGCVol', 'V1surface', 'V1thickness'};
names = {'Optic Tract FC','Optic Tract FD','fdc', 'FA', 'MD', 'Optic Radiation FC', 'Optic Radiation FD', 'Optic Radiation FDC', 'Intracranial Volume', 'LGN Volume', 'Mean Adjusted RGC Volume', 'V1 Surface', 'V1 Thickness'};
measureSet = {'Height_inches','Weight_pounds', 'Axial_Length_average', 'TIV'};

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
        
        % Make bar plots for correlations
        if strcmp(measureSet{ii}, 'Height_inches')
            ValheightR = R(1,2);
            ValheightConfidence = std(bootR);
        end
        if strcmp(measureSet{ii}, 'Weight_pounds')
            ValweightR = R(1,2);
            ValweightConfidence = std(bootR);  
        end
        if strcmp(measureSet{ii}, 'TIV')
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
    ylim([-1 1])
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