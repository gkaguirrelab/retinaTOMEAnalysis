function comboTable=lineThicknessAnalysis(GCIPthicknessFile, varargin)
% Do some analysis
%
% Description:
%   Foo
%
% Examples:
%{
    dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
    GCIPthicknessFile = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg');
    comboTable=lineThicknessAnalysis(GCIPthicknessFile);
%}

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('GCIPthicknessFile',@ischar);

% Optional analysis params
p.addParameter('showPlots',true,@islogical);
p.addParameter('dataSaveName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'), 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
p.addParameter('mmPerDegFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat'),@ischar);
p.addParameter('figSaveDir','/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/horizontalLine',@ischar);
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);



%% Parse and check the parameters
p.parse(GCIPthicknessFile, varargin{:});

loadAndMergeData


convertThicknessToVolume

%% Save Directory

saveDir = 'D:\Min\Dropbox (Personal)\Research\Projects\retinaTOMEAnalysis\code\figures\PaperFigures2';
mkdir(saveDir);
for i=1:9
    mkdir(fullfile(saveDir,['fig' num2str(i)]));
end


%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVec)'),num2cell(nanmean(gcVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcMeanThick','gcVolumePerDegSq'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');

fig1_ShowMontage
fig2_ThickAndVolRelationships

%% Conduct PCA upon the tissue volume data

createVolumePCA

% Show the effect of smoothing on the PCA scores
nDimsToUse = 6;

fig3_EffectOfSmoothingOnPCA


%% Identify the scores that correlate with axial length
%%%%%%%%%%%%%%%
%%
%% NEED TO ORTHOGONALIZE X WITH RESPECT TO AXIAL LENGTH HERE
%%
%% AND ADD BACK IN THE EFFECT OF THE EMMETROPIC EYE (23.58 mm)
%%
%% also synthesize
adjustAndSynthPCAWithAxialLength

fig4_smoothedPCAReconstruction

fig5_PCACompRegressionAxialLength

%%%%%%%%%%%%%%%
% Plot the reconstructions with the adjustment
fig6_AxialLengthAdjustedReconstruction

%%%%%%%%%%%%%%%
% Plot the synthesized reconstructions by axial length
fig7_SynthesizeALImpactOnEmmEye

% %%%%%%%%%%%%%%%
%figure 8 profiles after reconstruction with the AL component removed
fig8_VolRelationshipAfterAdjustment


%figure 9 profiles after reconstruction using the synthetic AL component
fig9_SynthProfileAndALRelationship

%% Save some results
%save(p.Results.dataSaveName,'XPos_Degs','meanPC1','gcVolumePerDegSqAdjust','coeff','scoreExpandedSmoothed');

end

