
% Obtain the dropBox base directory
dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');



saveDir = '/Users/aguirre/Documents/MATLAB/projects/rgcPopulationModel/data/';
%
%input directories
inputVolandSegDir = fullfile(dropboxBaseDir,'AOSO_analysis','OCTExplorerSegmentationData');

%output directories
thicknessMapSaveDir = fullfile(dropboxBaseDir,'AOSO_analysis','2DThicknessMapsAllLayers_MinChenMontage');
rgcMapSavePath = fullfile(saveDir,'rgcIplThicknessMap.mat');
resultTableFileName = fullfile(saveDir,'data/octRGCResultTable.csv');
 
mkdir(thicknessMapSaveDir);
mkdir(saveDir);

makeHorVertOCTMontage(inputVolandSegDir,thicknessMapSaveDir)
analyzeThicknessMaps(thicknessMapSaveDir, rgcMapSavePath, resultTableFileName)
