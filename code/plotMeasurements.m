%% plotMeasurements

% This script plots the measurements for a subject from the four meridians
% (nasal, temporal, superior, inferior) all on the same x-axis of
% eccentricity values.

% Obtain the dropBox base directory
dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');

%% load & place data

% Find and load the files we need
horizThicknessProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', ...
    'OCTExplorerExtendedHorizontalData', 'GCIP_thicknessesByDeg.mat');
load(horizThicknessProfile,'GCthicknessValuesAtXPos_um','XPos_Degs');

% Place data into a structure

thickData(1).label = 'temporal';
thickData(1).angle = 180;
thickData(1).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
thickData(1).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs<=0)')';

thickData(2).label = 'nasal';
thickData(2).angle = 0;
thickData(2).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
thickData(2).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs>=0)')';

% Find and load the files we need
vertThicknessProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', ...
    'OCTSingleVerticalData', 'GCIP_thicknessesByDeg.mat');
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um','XPos_Degs');

% Place data into a structure

thickData(3).label = 'superior';
thickData(3).angle = 90;
thickData(3).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
thickData(3).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs<=0)')';

thickData(4).label = 'inferior';
thickData(4).angle = 270;
thickData(4).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
thickData(4).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs>=0)')';

%% plot measurements
plot(XPos_Degs, squeeze(GCthicknessValuesAtXPos_um(1,1,:)))
