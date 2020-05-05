%%%
%Generates figures for Papers
close all
%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('GCIPthicknessFile',@ischar);

% Optional analysis params
p.addParameter('showPlots',true,@islogical);
%dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
% p.addParameter('dataSaveName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'), 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
% p.addParameter('mmPerDegFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat'),@ischar);
% p.addParameter('figSaveDir','/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/horizontalLine',@ischar);
% p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

dropboxBaseDir ='D:\Min\Dropbox (Aguirre-Brainard Lab)';
p.addParameter('dataSaveName',fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
p.addParameter('mmPerDegFileName',fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat'),@ischar);
p.addParameter('subjectTableFileName',fullfile(dropboxBaseDir,'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

saveDir = 'D:\Min\Dropbox (Personal)\Research\Projects\retinaTOMEAnalysis\code\figures\PaperFigures2';
mkdir(saveDir);
for i=1:9
    mkdir(fullfile(saveDir,['fig' num2str(i)]));
end
%% Parse and check the parameters
GCIPthicknessFile = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg.mat');
p.parse(GCIPthicknessFile);

loadAndMergeData

%% Convert from mm thickness to tissue volume

convertThicknessToVolume


%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVec)'),num2cell(nanmean(gcVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcMeanThick','gcVolumePerDegSq'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');


%% Conduct PCA upon the tissue volume data

createVolumePCA



%% Identify the scores that correlate with axial length


%%%%%%%%%%%%%%%
%%
%% ORTHOGONALIZE X WITH RESPECT TO AXIAL LENGTH
%% AND ADD BACK IN THE EFFECT OF THE EMMETROPIC EYE (23.58 mm)
%% AND synthesize
adjustAndSynthPCAWithAxialLength

%Figure 1 - Shows the montage and associated parts.
fig1_ShowMontage

%Figure 2
fig2_ThickAndVolRelationships

%Figure 3
fig3_EffectOfSmoothingOnPCA

%Figure 4
fig4_smoothedPCAReconstruction

%Figure 5
fig5_PCACompRegressionAxialLength

%Figure 6
%%%%%%%%%%%%%%%

fig6_AxialLengthAdjustedReconstruction

%Figure 7 
fig7_SynthesizeALImpactOnEmmEye


% %%%%%%%%%%%%%%%
%figure 8 profiles after reconstruction with the AL component removed
fig8_VolRelationshipAfterAdjustment


%figure 9 profiles after reconstruction using the synthetic AL component
fig9_SynthProfileAndALRelationship

% figure
% profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(1,1:nDimsToUse)';
% plot(profileFit,'-r','LineWidth',1);
% hold on
% profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(16,1:nDimsToUse)';
% plot(profileFit,'-g','LineWidth',1);
% profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(50,1:nDimsToUse)';
% plot(profileFit,'-b','LineWidth',1);
% ylim([-1 8])
% title('Synthetic from Axial Lengths: Red=21.79, Green=23.58, Blue=27.57')






%% Save some results
%save(p.Results.dataSaveName,'XPos_Degs','meanPC1','gcVolumePerDegSqAdjust','coeff','scoreExpandedSmoothed');


