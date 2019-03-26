addpath('.\SupportFunctions');

%baseDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/';
baseDir = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)'
saveDir = '/Users/aguirre/Documents/MATLAB/projects/rgcPopulationModel/data/';
saveDir = 'C:\Users\dontm\Downloads\temp';
%
%input directories
inputVolandSegDir = fullfile(baseDir,'AOSO_analysis','OCTExplorerSegmentationData');

%output directories
thicknessMapSaveDir = fullfile(baseDir,'AOSO_analysis','2DThicknessMapsAllLayers_MinChenMontage');
thicknessMapSaveDir = 'C:\Users\dontm\Downloads\temp';
rgcMapSavePath = fullfile(saveDir,'rgcIplThicknessMap.mat');
resultTableFileName = fullfile(saveDir,'data/octRGCResultTable.csv');
 
mkdir(thicknessMapSaveDir);
mkdir(saveDir);

generateAlignedThicknessMaps(inputVolandSegDir,thicknessMapSaveDir)
analyzeThicknessMaps(thicknessMapSaveDir, rgcMapSavePath, resultTableFileName)
