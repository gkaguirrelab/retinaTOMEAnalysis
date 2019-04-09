
% Obtain the dropBox base directory
dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');


%% octExplorerXMLToMatlab
% The OCT data are exported from the Spectralis system in '.vol' format,
% and then processed using OCT Explorer v5.0 (on a Macintosh). This
% analysis yields (for each '.vol' file) a file named
% `_Surfaces_Retina-JEI-Final.xml`. This first step converts these XML
% files to a format more easily loaded in matlab.
dataDir = fullfile(dropboxBaseDir,'AOSO_analysis','OCTExplorerSegmentationData');
octExplorerXMLToMatlab(dataDir);


%% makeHorizontalVerticalMontage
% Each eye was studied with a macular volume scan that was extended in the
% horizontal and in the vertical direction. This step creates an aligned
% montage of these two measurements.
dataDir = fullfile(dropboxBaseDir,'AOSO_analysis','OCTExplorerSegmentationData');
saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','2DThicknessMapsAllLayers_MinChenMontage');
makeHorizontalVerticalMontage(dataDir,saveDir)

%% makeleftRightMontage

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
