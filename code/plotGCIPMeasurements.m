%% plotGCIPMeasurements

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
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');

subList = {};
thickVec = [];
gcVec = [];
ipVec = [];
ratioVec = [];
gcipMeanOD = [];
gcipMeanOS = [];

% Obtain the GC thickness and ratio functions for each subject. While we
% are at it, confirm that there is a substantial correlation across
% subjects between the left and right eye in the median of the ratio
% functions.
GCthicknessValuesAtXPos_um(GCthicknessValuesAtXPos_um==0)=nan;
IPthicknessValuesAtXPos_um(IPthicknessValuesAtXPos_um==0)=nan;

for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % We are keeping this subject
        subList(end+1) = {subIDs(ii,:)};
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,:))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,:)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,:))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,:)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec(:,end+1) = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec(:,end+1) = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            fprintf([subList{end} ': missing eye\n']);
            gcVec(:,end+1) = nanmean([gcVecOD,gcVecOS],2);
            ipVec(:,end+1) = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the ratio and thickness vecs
        thickVec(:,end+1) = sum([gcVec(:,end),ipVec(:,end)],2,'includenan');
        ratioVec(:,end+1) = gcVec(:,end)./thickVec(:,end);
        
        % Save the m value for each eye and layer
        gcipMeanOD(end+1) = nanmean(gcVecOD+ipVecOD);
        gcipMeanOS(end+1) = nanmean(gcVecOS+ipVecOS);
        
    end
end

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanThickVecProfile = nanmean(thickVec,2);
meanGCVecProfile = nanmean(gcVec,2);
meanIPVecProfile = nanmean(ipVec,2);
meanRatioVecProfile = nanmean(ratioVec,2);

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanGCVecProfile(badIdx)=nan;

% plot horizontal GC thickness
subplot(2,2,1)
plot(XPos_Degs, squeeze(meanGCVecProfile))
hold on 

% plot horizontal IP thickness
plot(XPos_Degs, squeeze(meanIPVecProfile))
title('Average Thickness Along the Horizontal Meridian')
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
hold off

% plot horizontal GC:GCIP thickness
subplot(2,2,3)
plot(XPos_Degs, squeeze(meanRatioVecProfile))
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

subList = {};
thickVec = [];
gcVec = [];
ipVec = [];
ratioVec = [];
gcipMeanOD = [];
gcipMeanOS = [];

% Obtain the GC thickness and ratio functions for each subject. While we
% are at it, confirm that there is a substantial correlation across
% subjects between the left and right eye in the median of the ratio
% functions.
GCthicknessValuesAtXPos_um(GCthicknessValuesAtXPos_um==0)=nan;
IPthicknessValuesAtXPos_um(IPthicknessValuesAtXPos_um==0)=nan;

for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % We are keeping this subject
        subList(end+1) = {subIDs(ii,:)};
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,:))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,:)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,:))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,:)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec(:,end+1) = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec(:,end+1) = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            fprintf([subList{end} ': missing eye\n']);
            gcVec(:,end+1) = nanmean([gcVecOD,gcVecOS],2);
            ipVec(:,end+1) = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the ratio and thickness vecs
        thickVec(:,end+1) = sum([gcVec(:,end),ipVec(:,end)],2,'includenan');
        ratioVec(:,end+1) = gcVec(:,end)./thickVec(:,end);
        
        % Save the m value for each eye and layer
        gcipMeanOD(end+1) = nanmean(gcVecOD+ipVecOD);
        gcipMeanOS(end+1) = nanmean(gcVecOS+ipVecOS);
        
    end
end

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanThickVecProfile = nanmean(thickVec,2);
meanGCVecProfile = nanmean(gcVec,2);
meanIPVecProfile = nanmean(ipVec,2);
meanRatioVecProfile = nanmean(ratioVec,2);

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanGCVecProfile(badIdx)=nan;

% plot vertical GC thickness
subplot(2,2,2)
plot(XPos_Degs, squeeze(meanGCVecProfile))
hold on

% plot vertical IP thickness
plot(XPos_Degs, squeeze(meanIPVecProfile))
title('Average Thickness Along the Vertical Meridian')
legend('GCL','IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
hold off

% plot horizontal GC:GCIP thickness
subplot(2,2,4)
plot(XPos_Degs, squeeze(meanRatioVecProfile))
title('Average GCL Ratio Along the Vertical Meridian')
legend('GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])