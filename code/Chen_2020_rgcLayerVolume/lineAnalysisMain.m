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

dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));
%dropboxBaseDir ='C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)';



%% Create Save Directory
saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','GCPaperFigures');
mkdir(saveDir);
for i=1:9
    mkdir(fullfile(saveDir,['fig' num2str(i)]));
end


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('GCIPthicknessFile',@ischar);

% Optional analysis params
p.addParameter('showPlots',true,@islogical);
p.addParameter('dataSaveName',fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
p.addParameter('mmPerDegFileName',fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat'),@ischar);
p.addParameter('figSaveDir',saveDir,@ischar);
p.addParameter('subjectTableFileName',fullfile(dropboxBaseDir,'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);



%% Parse and check the parameters
p.parse(GCIPthicknessFile, varargin{:});

%loads the GC profiles, mean gc profile, bad indexes, subject list, and X positions
[gcVec,meanGCVecProfile,badIdx,subList,XPos_Degs,subjectTable,thicknessTable] =  loadAndMergeData(p,GCIPthicknessFile);

%Converts GC thickness into GC volume and provide mmSqPerDegSq conversion for each subject
[mmSqPerDegSq,gcVolumePerDegSq,meanGCVolumePerDegSqProfile,volumeTable] = convertThicknessToVolume(p,gcVec,badIdx,subList,XPos_Degs,subjectTable);

% Join the thickness and volume data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
gcTable = join(thicknessTable,volumeTable,'Keys','AOSO_ID');
comboTable = join(gcTable,subjectTable,'Keys','AOSO_ID');


nDimsToUse = 6;%number of PCA components to use.

%% Conduct PCA upon the tissue volume data
%Perform the PCA and smooth components
[GCVolPCAScoreExpanded, GCVolPCAScoreExpandedSmoothed, GCVolPCACoeff, GCVolPCAVarExplained] = createVolumePCA(gcVolumePerDegSq,badIdx,XPos_Degs);
%Adjust each coeff with the Axial length contribution, also create some synthetic ones
[adjustedGCVolPCACoeff,synGCVolPCACoeff,synthALRange] = adjustAndSynthPCAWithAxialLength(nDimsToUse,GCVolPCACoeff,comboTable.Axial_Length_average);


if(p.Results.showPlots)%Flag determines if we plot and save or not
    close all;
    
    %Pulls up an example montage and individual pieces
    fig1_ShowMontage(dropboxBaseDir, saveDir)
    
    %Plots Relating GC thickness and volume to axial length
    fig2_ThickAndVolRelationships(XPos_Degs, gcVec, meanGCVecProfile, mmSqPerDegSq, gcVolumePerDegSq, meanGCVolumePerDegSqProfile,comboTable,saveDir)
    
    % Panel d of figure 2 is a couple of extreme model eyes
    fig2d_ExtremeEyes(saveDir)
    
    %Show the effect of smoothing on each of the PCA scores
    fig3_EffectOfSmoothingOnPCA(GCVolPCAVarExplained,GCVolPCAScoreExpanded,GCVolPCAScoreExpandedSmoothed,nDimsToUse,saveDir)
    
    %reconstruct original profiles using smoothed and first 'nDimsToUse'
    %components to show they still looks correct
    fig4_smoothedPCAReconstruction(gcVolumePerDegSq,GCVolPCAScoreExpandedSmoothed,GCVolPCACoeff,nDimsToUse,saveDir)
    
    %Look at how each PCA component regresses with axial length
    fig5_PCACompRegressionAxialLength(comboTable.Axial_Length_average,GCVolPCACoeff,nDimsToUse,saveDir)
    
    % Plot the reconstructions of each of the profiles after adjusting the PCs
    % to remove the influence of axial length (subjects sorted by axial length)
    fig6_AxialLengthAdjustedReconstruction(comboTable.Axial_Length_average,gcVolumePerDegSq,GCVolPCAScoreExpandedSmoothed,adjustedGCVolPCACoeff,nDimsToUse,saveDir)
    
    % Using PCa, show how axial length will affect the average emmetropic profile
    fig7_SynthesizeALImpactOnEmmEye(GCVolPCAScoreExpandedSmoothed,synGCVolPCACoeff,XPos_Degs,saveDir)
    
    %Show AL adjusted profiles stacked, and their mean/median relationship with axial length before/after adjustment
    fig8_VolRelationshipAfterAdjustment(XPos_Degs,comboTable,GCVolPCAScoreExpandedSmoothed,adjustedGCVolPCACoeff,GCVolPCACoeff,nDimsToUse,saveDir)
    
    %Synthetic profiles demonstrating only the influence of the axial length
    %component i PC space
    fig9_SynthProfileAndALRelationship(XPos_Degs,GCVolPCAScoreExpandedSmoothed,synGCVolPCACoeff,[min(synthALRange) 23.58 max(synthALRange)],nDimsToUse,saveDir)
end

end

