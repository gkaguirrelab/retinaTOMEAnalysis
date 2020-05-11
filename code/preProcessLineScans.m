%% preProcessLineScans
%
% This script performs initial processing on the vertical and horizontal
% line scans collected using OCT as part of the TOME project.
%

% Obtain the dropBox base directory
dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');


%% Convert vertical line scans
% We collected a single, vertical line scan through the fovea. These images
% were exported from the Heidelberg Spectralis device in ".vol" format.
% This first step converts the files to ".nii" format, and saves these in
% the analysis directory.
inputDir = fullfile(dropboxBaseDir,'AOSO_data','connectomeRetinaData');
outputDir =fullfile(dropboxBaseDir, 'AOSO_analysis','OCTSingleVerticalData');
saveSingleVerticalOCT(inputDir,outputDir)


%% Convert and montage horizontal line scans
% We collected three, overapping horizontal line scans. Each line scan
% included a portion of the fovea. At this stage the three line scans are
% montaged into a single representation, and saved in ".nii" format.
inputDir = fullfile(dropboxBaseDir,'AOSO_data','connectomeRetinaData');
outputDir =fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');

% The scanInfoFile is an Excel file that has information regarding which
% OCT image is the left, right, or center in the montage. This information
% is needed to guide the montage process.
scanInfoFile =fullfile(outputDir,'ExtendedOCTDataDescriptions.xlsx');

% Peform the analysis
montageExtendedOCT(inputDir,scanInfoFile,outputDir)


%% Smooth the manual segmentations 
% After the line scans are converted to nii format, a human operator
% uses ITKSNAP to perform manual labeling of the GC and IP layers. This
% next stage loads the labels and fits smooth splines to the three borders
% implied by the two labels.
inputDir = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');
smoothManualSegGCIPExtOCT(inputDir)


%% Caclulate IP and GC Thickness Values
% This stage calculates and saves a .mat file with the thickness of each of
% the two layers defined by the three boundaries.
inputDir = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData');
calcExtOCTGCIPThickness(inputDir)


%% Create the mmSquarePerDegreeSquare maps
% We create for each subject a model eye based upon their biometric
% measurements. We calculate for each model eye the square mm of retina
% that subtend a square degree of visual angle at several points. This set
% of measurements is then fit with a polynomial surface. The polynomial
% allows us to convert mm2 to degree2 for any arbitrary position in the
% retina. The set of polynomial fits are stored in a single matlab variable
% called "mmPerDegPolyFit.mat".
outputDir = fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps');
makeMmPerDegMaps(saveDir)


%% Analyze GC IP layers in the extended OCT 
GCIPthicknessFile = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg.mat');
comboTable=lineThicknessAnalysis(GCIPthicknessFile);