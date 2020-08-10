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
horizRatioVec = [];
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
        horizRatioVec(:,end+1) = gcVec(:,end)./thickVec(:,end);
        
        % Save the m value for each eye and layer
        gcipMeanOD(end+1) = nanmean(gcVecOD+ipVecOD);
        gcipMeanOS(end+1) = nanmean(gcVecOS+ipVecOS);
        
    end        
end

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanHorizThickVecProfile = nanmean(thickVec,2);
meanHorizGCVecProfile = nanmean(gcVec,2);
meanHorizIPVecProfile = nanmean(ipVec,2);
meanHorizRatioVecProfile = nanmean(horizRatioVec,2);

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanHorizGCVecProfile(badIdx)=nan;
meanHorizIPVecProfile(badIdx)=nan;
meanHorizRatioVecProfile(badIdx)=nan;
meanHorizThickVecProfile(badIdx)=nan;

% plot horizontal GC thickness
subplot(2,2,1)
plot(XPos_Degs, squeeze(meanHorizGCVecProfile))
hold on

% plot outliers 11069 and 11076
if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(22,1,:)))) || ...
        ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(22,2,:))))

    % Get the data for each layer and eye and convert to mm
    gcVecOD = squeeze(GCthicknessValuesAtXPos_um(22,1,:))/1000;
    gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(22,2,:)))/1000;
    ipVecOD = squeeze(IPthicknessValuesAtXPos_um(22,1,:))/1000;
    ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(22,2,:)))/1000;

    % Detect if the data from one eye is missing
    if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
        subOneGcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        subOneIpVec = mean([ipVecOD,ipVecOS],2,'includenan');
    else
        subOneGcVec = nanmean([gcVecOD,gcVecOS],2);
        subOneIpVec = nanmean([ipVecOD,ipVecOS],2);
    end

    % Calculate the thickness vec
    subOneThickVec = sum([subOneGcVec,subOneIpVec],2,'includenan');

    % plot subject
    plot(XPos_Degs, subOneGcVec);
    plot(XPos_Degs, subOneIpVec);
end

if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(29,1,:)))) || ...
        ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(29,2,:))))

    % Get the data for each layer and eye and convert to mm
    gcVecOD = squeeze(GCthicknessValuesAtXPos_um(29,1,:))/1000;
    gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(29,2,:)))/1000;
    ipVecOD = squeeze(IPthicknessValuesAtXPos_um(29,1,:))/1000;
    ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(29,2,:)))/1000;

    % Detect if the data from one eye is missing
    if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
        subTwoGcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        subTwoIpVec = mean([ipVecOD,ipVecOS],2,'includenan');
    else
        subTwoGcVec = nanmean([gcVecOD,gcVecOS],2);
        subTwoIpVec = nanmean([ipVecOD,ipVecOS],2);
    end

    % Calculate the thickness vec
    subTwoThickVec = sum([subTwoGcVec,subTwoIpVec],2,'includenan');

    % plot subject
    plot(XPos_Degs, subTwoGcVec);
    plot(XPos_Degs, subTwoIpVec);
end

% plot horizontal IP thickness
plot(XPos_Degs, squeeze(meanHorizIPVecProfile))
title('Average Thickness Along the Horizontal Meridian')
legend('GCL', '11069 GCL', '11069 IPL', '11076 GCL', '11076 IPL', 'IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
hold off

% plot horizontal GC:GCIP thickness
subplot(2,2,3)
plot(XPos_Degs, squeeze(meanHorizRatioVecProfile))
hold on
plot(XPos_Degs, squeeze(subOneGcVec ./ subOneThickVec))
plot(XPos_Degs, squeeze(subTwoGcVec ./ subTwoThickVec))
title('Average GCL Ratio Along the Horizontal Meridian')
legend('GC:GCIP Average', '11069 GC:GCIP', '11076 GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])
hold off


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
vertRatioVec = [];
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
        vertRatioVec(:,end+1) = gcVec(:,end)./thickVec(:,end);
        
        % Save the m value for each eye and layer
        gcipMeanOD(end+1) = nanmean(gcVecOD+ipVecOD);
        gcipMeanOS(end+1) = nanmean(gcVecOS+ipVecOS);
        
    end
end

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanVertThickVecProfile = nanmean(thickVec,2);
meanVertGCVecProfile = nanmean(gcVec,2);
meanVertIPVecProfile = nanmean(ipVec,2);
meanVertRatioVecProfile = nanmean(vertRatioVec,2);

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanVertGCVecProfile(badIdx)=nan;
meanVertIPVecProfile(badIdx)=nan;
meanVertRatioVecProfile(badIdx)=nan;
meanVertThickVecProfile(badIdx)=nan;

% plot vertical GC thickness
subplot(2,2,2)
plot(XPos_Degs, squeeze(meanVertGCVecProfile))
hold on

% plot outliers 11069 and 11076
if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(22,1,:)))) || ...
        ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(22,2,:))))

    % Get the data for each layer and eye and convert to mm
    gcVecOD = squeeze(GCthicknessValuesAtXPos_um(22,1,:))/1000;
    gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(22,2,:)))/1000;
    ipVecOD = squeeze(IPthicknessValuesAtXPos_um(22,1,:))/1000;
    ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(22,2,:)))/1000;

    % Detect if the data from one eye is missing
    if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
        subOneGcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        subOneIpVec = mean([ipVecOD,ipVecOS],2,'includenan');
    else
        subOneGcVec = nanmean([gcVecOD,gcVecOS],2);
        subOneIpVec = nanmean([ipVecOD,ipVecOS],2);
    end

    % Calculate the thickness vec
    subOneThickVec = sum([subOneGcVec,subOneIpVec],2,'includenan');

    % plot subject
    plot(XPos_Degs, subOneGcVec);
    plot(XPos_Degs, subOneIpVec);
end

if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(29,1,:)))) || ...
        ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(29,2,:))))

    % Get the data for each layer and eye and convert to mm
    gcVecOD = squeeze(GCthicknessValuesAtXPos_um(29,1,:))/1000;
    gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(29,2,:)))/1000;
    ipVecOD = squeeze(IPthicknessValuesAtXPos_um(29,1,:))/1000;
    ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(29,2,:)))/1000;

    % Detect if the data from one eye is missing
    if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
        subTwoGcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        subTwoIpVec = mean([ipVecOD,ipVecOS],2,'includenan');
    else
        subTwoGcVec = nanmean([gcVecOD,gcVecOS],2);
        subTwoIpVec = nanmean([ipVecOD,ipVecOS],2);
    end

    % Calculate the thickness vec
    subTwoThickVec = sum([subTwoGcVec,subTwoIpVec],2,'includenan');

    % plot subject
    plot(XPos_Degs, subTwoGcVec);
    plot(XPos_Degs, subTwoIpVec);
end

% plot vertical IP thickness
plot(XPos_Degs, squeeze(meanVertIPVecProfile))
title('Average Thickness Along the Vertical Meridian')
legend('GCL', '11069 GCL', '11069 IPL', '11076 GCL', '11076 IPL', 'IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
hold off

% plot vertical GC:GCIP thickness
subplot(2,2,4)
plot(XPos_Degs, squeeze(meanVertRatioVecProfile))
hold on
plot(XPos_Degs, squeeze(subOneGcVec ./ subOneThickVec))
plot(XPos_Degs, squeeze(subTwoGcVec ./ subTwoThickVec))
title('Average GCL Ratio Along the Vertical Meridian')
legend('GC:GCIP Average', '11069 GC:GCIP', '11076 GC:GCIP')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])
hold off

%% plot GC:GCIP across all four meridians
% % this section has a bug in it in that the meridians are flipped in left
% % vs. right eyes
% f2 = figure;
% 
% % create variables for all four meridans (leaves out value at idx 1281 to
% % keep arrays even, which is at a eccentricity of 0 so there is no GCIP there)
% meanTempGCVecProfile = flipud(meanHorizGCVecProfile(1:1280));
% meanTempThickVecProfile = flipud(meanHorizThickVecProfile(1:1280));
% meanAbsTempThickVecProfile = meanTempGCVecProfile ./ (meanTempThickVecProfile .^ 2);
% meanNasalGCVecProfile = meanHorizGCVecProfile(1282:2561);
% meanNasalThickVecProfile = meanHorizThickVecProfile(1282:2561);
% meanAbsNasalThickVecProfile = meanNasalGCVecProfile ./ (meanNasalThickVecProfile .^ 2);
% meanSupGCVecProfile = flipud(meanVertGCVecProfile(1:1280));
% meanSupThickVecProfile = flipud(meanVertThickVecProfile(1:1280));
% meanAbsSupThickVecProfile = meanSupGCVecProfile ./ (meanSupThickVecProfile .^ 2);
% meanInfGCVecProfile = meanVertGCVecProfile(1282:2561);
% meanInfThickVecProfile = meanVertThickVecProfile(1282:2561);
% meanAbsInfThickVecProfile = meanInfGCVecProfile ./ (meanInfThickVecProfile .^ 2);
% 
% % get averages across the four meridians
% meanGCVecProfile = meanTempGCVecProfile + meanNasalGCVecProfile + ...
%     meanSupGCVecProfile + meanInfGCVecProfile;
% meanGCVecProfile = meanGCVecProfile ./ 4;
% meanThickVecProfile = meanTempThickVecProfile + meanNasalThickVecProfile + ...
%     meanSupThickVecProfile + meanInfThickVecProfile;
% meanThickVecProfile = meanThickVecProfile ./ 4;
% meanRatioVecProfile = meanGCVecProfile ./ meanThickVecProfile;
% 
% posEccentricity = XPos_Degs(1282:2561);
% 
% % plot
% subplot(1, 2, 1)
% plot(posEccentricity, squeeze(meanRatioVecProfile))
% title('Average GCL Ratio Along the Four Meridians')
% legend('GC:GCIP')
% xlabel('Eccentricity') 
% ylabel('GCL:GCIPL Ratio')
% xlim([0 30])
% ylim([0 .75])
% 
% subplot(1, 2, 2)
% plot(posEccentricity, squeeze(meanAbsTempThickVecProfile))
% hold on
% plot(posEccentricity, squeeze(meanAbsNasalThickVecProfile))
% plot(posEccentricity, squeeze(meanAbsSupThickVecProfile))
% plot(posEccentricity, squeeze(meanAbsInfThickVecProfile))
% title('Average GCL Ratio Along the Four Meridians')
% legend('Temporal', 'Nasal', 'Superior', 'Inferior')
% xlabel('Eccentricity') 
% ylabel('GCL:GCIPL^2 Ratio')
% xlim([0 30])
% hold off

%% plot individual subjects
f3 = figure;

load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');
posEccentricity = XPos_Degs(1282:2561);

subplot(2, 2, 1)
% plot temporal meridian
for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,1:1280))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,1282:2561)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,1:1280))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,1282:2561)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            gcVec = nanmean([gcVecOD,gcVecOS],2);
            ipVec = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the thickness vec
        thickVec = sum([gcVec,ipVec],2,'includenan');
        
        % plot subject
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);
        plot(posEccentricity, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Temporal Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])
hold off

% plot nasal meridian
subplot(2, 2, 2)

for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,1282:2561))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,1:1280)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,1282:2561))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,1:1280)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            gcVec = nanmean([gcVecOD,gcVecOS],2);
            ipVec = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the thickness vec
        thickVec = sum([gcVec,ipVec],2,'includenan');
        
        % plot subject
        plot(posEccentricity, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Nasal Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])
hold off

% plot superior meridian
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');
posEccentricity = XPos_Degs(1282:2561);

subplot(2, 2, 3)
for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,1:1280))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,1282:2561)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,1:1280))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,1282:2561)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            gcVec = nanmean([gcVecOD,gcVecOS],2);
            ipVec = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the thickness vec
        thickVec = sum([gcVec,ipVec],2,'includenan');
        
        % plot subject
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);
        plot(posEccentricity, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Superior Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])
hold off

% plot inferior meridian
subplot(2, 2, 4)

for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,1,:)))) || ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(ii,2,:))))
        
        % Get the data for each layer and eye and convert to mm
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,1282:2561))/1000;
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,1:1280)))/1000;
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,1282:2561))/1000;
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,1:1280)))/1000;
        
        % Detect if the data from one eye is missing
        if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
            gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
            ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
        else
            gcVec = nanmean([gcVecOD,gcVecOS],2);
            ipVec = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the thickness vec
        thickVec = sum([gcVec,ipVec],2,'includenan');
        
        % plot subject
        plot(posEccentricity, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Inferior Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])
hold off

%% plot single meridian for an individual subject
f4 = figure;

load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');
posEccentricity = XPos_Degs(1282:2561);

subject = 2;

% plot temporal meridian
if ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(subject,1,:)))) || ...
        ~all(isnan(squeeze(GCthicknessValuesAtXPos_um(subject,2,:))))

    % Get the data for each layer and eye and convert to mm
    gcVecOD = squeeze(GCthicknessValuesAtXPos_um(subject,1,1:1280))/1000;
    gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(subject,2,1282:2561)))/1000;
    ipVecOD = squeeze(IPthicknessValuesAtXPos_um(subject,1,1:1280))/1000;
    ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(subject,2,1282:2561)))/1000;

    % Detect if the data from one eye is missing
    if ~all(isnan(gcVecOD)) && ~all(isnan(gcVecOS))
        gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
    else
        gcVec = nanmean([gcVecOD,gcVecOS],2);
        ipVec = nanmean([ipVecOD,ipVecOS],2);
    end

    % Calculate the thickness vec
    thickVec = sum([gcVec,ipVec],2,'includenan');
    
    gcVec = flipud(gcVec);
    thickVec = flipud(thickVec);
    ipVec = flipud(ipVec);

    % plot subject
    subplot(2, 2, 1)
    plot(posEccentricity, gcVec);
    title('GCL Thickness')
    xlabel('Eccentricity') 
    ylabel('Thickness')
    ylim([0 .07])
    subplot(2, 2, 2)
    plot(posEccentricity, ipVec);
    title('IPL Thickness')
    xlabel('Eccentricity') 
    ylabel('Thickness')
    ylim([0 .07])
    subplot(2, 2, 3)
    plot(posEccentricity, gcVec ./ thickVec);
    title('GC:GCIP Thickness')
    xlabel('Eccentricity') 
    ylabel('Thickness')
    subplot(2, 2, 4)
    plot(posEccentricity, gcVec ./ (thickVec .^2));
    title('GCL over GCIP')
    xlabel('Eccentricity') 
    ylabel('GC/(GCIP^2)')
end