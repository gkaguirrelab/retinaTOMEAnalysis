%% plotMeasurements

% This script plots the measurements for a subject from the four meridians
% (nasal, temporal, superior, inferior) all on the same x-axis of
% eccentricity values.

% Obtain the dropBox base directory
dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');   

% Instantiate a plotlab object
plotlabOBJ = plotlab();
    
% Apply the default plotlab recipe 
% overriding just the figure size
plotlabOBJ.applyRecipe(...
  'figureWidthInches', 20, ...
  'figureHeightInches', 15);

%% plot GCIP horizontal thickness

% Find and load the files we need across horizontal meridian
horizThicknessProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', ...
    'OCTExplorerExtendedHorizontalData', 'GCIP_thicknessesByDeg.mat');
load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');

% place horizontal GC thickness data into a structure
GCthickData(1).label = 'temporal';
GCthickData(1).angle = 180;
GCthickData(1).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
GCthickData(1).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs<=0)')';

GCthickData(2).label = 'nasal';
GCthickData(2).angle = 0;
GCthickData(2).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
GCthickData(2).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs>=0)')';

% plot horizontal GC thickness
subplot(2,2,1)
plot(XPos_Degs, squeeze(GCthicknessValuesAtXPos_um(1,1,:)))
hold on

% place horizontal IP thickness data into a structure
IPthickData(1).label = 'temporal';
IPthickData(1).angle = 180;
IPthickData(1).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
IPthickData(1).thickMM = (IPthicknessValuesAtXPos_um(XPos_Degs<=0)')';

IPthickData(2).label = 'nasal';
IPthickData(2).angle = 0;
IPthickData(2).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
IPthickData(2).thickMM = (IPthicknessValuesAtXPos_um(XPos_Degs>=0)')';

% plot horizontal IP thickness
plot(XPos_Degs, squeeze(IPthicknessValuesAtXPos_um(1,1,:)))
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness')
hold off

% plot horizontal GC:GCIP thickness
subplot(2,2,3)
plot(XPos_Degs, squeeze(GCthicknessValuesAtXPos_um(1,1,:)) ./ ... 
    (squeeze(IPthicknessValuesAtXPos_um(1,1,:)) + ... 
    squeeze(GCthicknessValuesAtXPos_um(1,1,:))))
legend('GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')


%% plot GCIP vertical thickness
% Find and load the files we need across the vertical meridian
vertThicknessProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', ...
    'OCTSingleVerticalData', 'GCIP_thicknessesByDeg.mat');
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');

% place vertical GC data into a structure
GCthickData(3).label = 'superior';
GCthickData(3).angle = 90;
GCthickData(3).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
GCthickData(3).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs<=0)')';

GCthickData(4).label = 'inferior';
GCthickData(4).angle = 270;
GCthickData(4).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
GCthickData(4).thickMM = (GCthicknessValuesAtXPos_um(XPos_Degs>=0)')';

% plot vertical GC thickness
subplot(2,2,2)
plot(XPos_Degs, squeeze(GCthicknessValuesAtXPos_um(1,1,:)))
hold on

% place vertical IP into a structure
IPthickData(3).label = 'superior';
IPthickData(3).angle = 90;
IPthickData(3).supportDeg = abs(XPos_Degs(XPos_Degs<=0));
IPthickData(3).thickMM = (IPthicknessValuesAtXPos_um(XPos_Degs<=0)')';

IPthickData(4).label = 'inferior';
IPthickData(4).angle = 270;
IPthickData(4).supportDeg = abs(XPos_Degs(XPos_Degs>=0));
IPthickData(4).thickMM = (IPthicknessValuesAtXPos_um(XPos_Degs>=0)')';

% plot vertical IP thickness
plot(XPos_Degs, squeeze(IPthicknessValuesAtXPos_um(1,1,:)))
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness')
hold off

% plot vertical GC:GCIP thickness
subplot(2,2,4)
plot(XPos_Degs, squeeze(GCthicknessValuesAtXPos_um(1,1,:)) ./ ... 
    (squeeze(IPthicknessValuesAtXPos_um(1,1,:)) + ... 
    squeeze(GCthicknessValuesAtXPos_um(1,1,:))))
legend('GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
