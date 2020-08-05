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

% plot horizontal GC thickness (flip the values in right eye so that the
% optic nerves for both eyes are located on the same x-values)
avgGC = mean(GCthicknessValuesAtXPos_um, 'omitnan');
leftEye = avgGC(1,1,:);
rightEye = flip(avgGC(1,2,:));
avgGC = (leftEye + rightEye) ./ 2;
subplot(2,2,1)
plot(XPos_Degs, squeeze(avgGC))
hold on

% plot horizontal IP thickness (flip the values in right eye so that the
% optic nerves for both eyes are located on the same x-values)
avgIP = mean(IPthicknessValuesAtXPos_um, 'omitnan');
leftEye = avgIP(1,1,:);
rightEye = flip(avgIP(1,2,:));
avgIP = (leftEye + rightEye) ./ 2;
plot(XPos_Degs, squeeze(avgIP))
title('Average Thickness Along the Horizontal Meridian')
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness')
xlim([-30 30])
ylim([0 70])
hold off

% plot horizontal GC:GCIP thickness (flip the values in right eye so that 
% the optic nerves for both eyes are located on the same x-values)
subplot(2,2,3)
avgGCIP = GCthicknessValuesAtXPos_um + IPthicknessValuesAtXPos_um;
avgGCIP = mean(avgGCIP, 'omitnan');
leftEye = avgGCIP(1,1,:);
rightEye = flip(avgGCIP(1,2,:));
avgGCIP = (leftEye + rightEye) ./ 2;
plot(XPos_Degs, squeeze(avgGC ./ avgGCIP))
title('Average GCL Ratio Along the Horizontal Meridian')
legend('GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])


%% plot GCIP vertical thickness

% Find and load the files we need across the vertical meridian
vertThicknessProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', ...
    'OCTSingleVerticalData', 'GCIP_thicknessesByDeg.mat');
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');

% plot vertical GC thickness (flip the values in right eye so that the
% optic nerves for both eyes are located on the same x-values)
avgGC = mean(GCthicknessValuesAtXPos_um, 'omitnan');
leftEye = avgGC(1,1,:);
rightEye = flip(avgGC(1,2,:));
avgGC = (leftEye + rightEye) ./ 2;
subplot(2,2,2)
plot(XPos_Degs, squeeze(avgGC))
hold on

% plot vertical IP thickness (flip the values in right eye so that the
% optic nerves for both eyes are located on the same x-values)
avgIP = mean(IPthicknessValuesAtXPos_um, 'omitnan');
leftEye = avgIP(1,1,:);
rightEye = flip(avgIP(1,2,:));
avgIP = (leftEye + rightEye) ./ 2;
plot(XPos_Degs, squeeze(avgIP))
plot(XPos_Degs, squeeze(avgIP))
title('Average Thickness Along the Vertical Meridian')
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness')
xlim([-30 30])
ylim([0 70])
hold off

% plot vertical GC:GCIP thickness (flip the values in right eye so that 
% the optic nerves for both eyes are located on the same x-values)
subplot(2,2,4)
avgGCIP = GCthicknessValuesAtXPos_um + IPthicknessValuesAtXPos_um;
avgGCIP = mean(avgGCIP, 'omitnan');
leftEye = avgGCIP(1,1,:);
rightEye = flip(avgGCIP(1,2,:));
avgGCIP = (leftEye + rightEye) ./ 2;
plot(XPos_Degs, squeeze(avgGC ./ avgGCIP))
title('Average GCL Ratio Along the Vertical Meridian')
legend('GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])
