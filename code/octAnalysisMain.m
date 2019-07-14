
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


%% makeLeftRightMontage
% The data from both eyes are positioned such that the fovea is at the
% center of the imaging space. The data from the left eye is mirror
% reversed. Data are averaged across the eyes. Signal values from close to
% the optic disc are trimmed away.
dataDir = fullfile(dropboxBaseDir,'AOSO_analysis','2DThicknessMapsAllLayers_MinChenMontage');
saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','averageThicknessMapsBySubject');
makeLeftRightMontage(dataDir,saveDir)



%% makeMmPerDegMaps
% Saves maps of the conversion of degree of visual angle to mm of retina
saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps');
makeMmPerDegMaps(saveDir)


%% makeVolumeMaps
% The maps of RGCIPL thickness are combined with the mmPerDeg maps to
% create maps of tissue volume per square visual degree.
thicknessMapDir = fullfile(dropboxBaseDir,'AOSO_analysis','averageThicknessMapsBySubject');
mmPerDegMapDir = fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps');
saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','volumeMapsBySubject');
makeVolumeMaps(thicknessMapDir,mmPerDegMapDir,saveDir);


%% Conduct a PCA analysis of the thickness maps
dataDir = fullfile(dropboxBaseDir,'AOSO_analysis','averageThicknessMapsBySubject');
thicknessMapPCAAnalysis(dataDir)


%% CODE TO PROCESS HORIZONTAL LINE SCANS

%% Montage Horizontal Extended OCT Scans
dataDir = fullfile(dropboxBaseDir,'AOSO_data','connectomeRetinaData');
saveDir =fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');
scanInfoFile =fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','ExtendedOCTDataDescriptions.xlsx');
montageExtendedOCT(dataDir,scanInfoFile,saveDir)

%% Smooth the manual segmentations 
dataDir = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');
smoothManualSegGCIPExtOCT(dataDir)

%% Analyze GC IP layers in the extended OCT 
dataDir = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');
analyzeGCIPExtendedOCT(dataDir)

