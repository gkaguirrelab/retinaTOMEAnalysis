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

subplot(2,2,3)
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
        
        % plot subject
        plot(XPos_Degs, gcVec(:,end) ./ thickVec(:,end), 'LineWidth', 1)
        hold on
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
horizbadIdx = subCountPerPoint<(length(subList)/3);
meanHorizGCVecProfile(horizbadIdx)=nan;
meanHorizIPVecProfile(horizbadIdx)=nan;
meanHorizRatioVecProfile(horizbadIdx)=nan;
meanHorizThickVecProfile(horizbadIdx)=nan;

% plot horizontal GC:GCIP thickness
plot(XPos_Degs, squeeze(meanHorizRatioVecProfile), 'color', 'black')
title('Average GCL Ratio Along the Horizontal Meridian')
% legend('GC:GCIP Average')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])
hold off

% plot horizontal GC thickness
subplot(2,2,1)
plot(XPos_Degs, squeeze(meanHorizGCVecProfile))
hold on

% plot horizontal IP thickness
plot(XPos_Degs, squeeze(meanHorizIPVecProfile))
title('Average Thickness Along the Horizontal Meridian')
legend('GCL', 'IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
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

subplot(2,2,4)
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
        
        % plot subject
        plot(XPos_Degs, gcVec(:,end) ./ thickVec(:,end), 'LineWidth', 1)
        hold on
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
vertbadIdx = subCountPerPoint<(length(subList)/3);
meanVertGCVecProfile(vertbadIdx)=nan;
meanVertIPVecProfile(vertbadIdx)=nan;
meanVertRatioVecProfile(vertbadIdx)=nan;
meanVertThickVecProfile(vertbadIdx)=nan;

% plot vertical GC:GCIP thickness
plot(XPos_Degs, squeeze(meanVertRatioVecProfile), 'color', 'black')
title('Average GCL Ratio Along the Vertical Meridian')
% legend('GC:GCIP Average')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL Ratio')
xlim([-30 30])
ylim([0 .75])

% plot vertical GC thickness
subplot(2,2,2)
plot(XPos_Degs, squeeze(meanVertGCVecProfile))
hold on

% plot vertical IP thickness
plot(XPos_Degs, squeeze(meanVertIPVecProfile))
title('Average Thickness Along the Vertical Meridian')
legend('GCL', 'IPL')
xlabel('Eccentricity') 
ylabel('Thickness (mm)')
xlim([-30 30])
ylim([0 .07])
hold off

%% plot individual subjects (distance)
f2 = figure;

load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');
posEccentricity = XPos_Degs(1282:2561);

mmProfile = fullfile(dropboxBaseDir, 'AOSO_analysis', 'mmProfile.mat');
load(mmProfile, 'XPos_mm');
xPos = XPos_mm{1};

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
        
        % get x-axis values
        xPosSub = flipud(xPos(1:1280, ii));
        xPosSub = abs(xPosSub);
        
        % plot subject
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);
        plot(xPosSub, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Temporal Meridian')
xlabel('Distance (mm)') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 10])
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
        xPosSub = xPos(1282:2561, ii);
        plot(xPosSub, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Nasal Meridian')
xlabel('Distance (mm)') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 10])
ylim([0 75])
hold off

% plot superior meridian
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');
posEccentricity = XPos_Degs(1282:2561);

xPos = XPos_mm{2};

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
        xPosSub = xPos(1282:2561, ii);
        plot(xPosSub, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Superior Meridian')
xlabel('Distance (mm)') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 5])
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
        
        % get x-axis values
        xPosSub = flipud(xPos(1:1280, ii));
        xPosSub = abs(xPosSub);
        
        % plot subject
        plot(xPosSub, gcVec ./ (thickVec .^2));
        hold on
    end
end
title('GCL Ratio Along the Inferior Meridian')
xlabel('Distance (mm)') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 5])
ylim([0 75])
hold off


%% plot standard deviation across subjects for each meridian
f3 = figure;

load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');
posEccentricity = XPos_Degs(1282:2561);

tempThickVec = [];
tempGcVec = [];
tempIpVec = [];
tempRatioVec = [];

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
        
        % Calculate the thickness vec and flip vecs
        thickVec = sum([gcVec,ipVec],2,'includenan');
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);
        
        tempGcVec(:, end+1) = gcVec;
        tempIpVec(:, end+1) = ipVec;
        tempThickVec(:, end+1) = thickVec;
        
        % Calculate the ratio and thickness vecs
        tempThickVec(:,end+1) = sum([tempGcVec(:,end),tempIpVec(:,end)],2,'includenan');
        tempRatioVec(:,end+1) = tempGcVec(:,end)./(tempThickVec(:,end) .^ 2);
    end
end
SD = std(tempRatioVec .', 'omitnan');
meanTempRatioVecProfile = nanmean(tempRatioVec,2);
errorbar(posEccentricity, meanTempRatioVecProfile, SD, '-s','MarkerSize',...
    2, 'MarkerEdgeColor','black','MarkerFaceColor','black')
% plot(posEccentricity, tempRatioVec)
title('Mean GC:GCIP^2 Along the Temporal Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])

% plot nasal meridian
nasThickVec = [];
nasGcVec = [];
nasIpVec = [];
nasRatioVec = [];

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
        
        nasGcVec(:, end+1) = gcVec;
        nasIpVec(:, end+1) = ipVec;
        nasThickVec(:, end+1) = thickVec;
        
        % Calculate the ratio and thickness vecs
        nasThickVec(:,end+1) = sum([nasGcVec(:,end),nasIpVec(:,end)],2,'includenan');
        nasRatioVec(:,end+1) = nasGcVec(:,end)./(nasThickVec(:,end) .^ 2);
    end
end
SD = std(nasRatioVec .', 'omitnan');
meanNasRatioVecProfile = nanmean(nasRatioVec,2);
errorbar(posEccentricity, meanNasRatioVecProfile, SD, '-s','MarkerSize',...
    2, 'MarkerEdgeColor','black','MarkerFaceColor','black')
title('Mean GC:GCIP^2 Along the Nasal Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 30])
ylim([0 75])

% plot superior meridian
load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs');
posEccentricity = XPos_Degs(1282:2561);

supThickVec = [];
supGcVec = [];
supIpVec = [];
supRatioVec = [];

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
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);
        
        supGcVec(:, end+1) = gcVec;
        supIpVec(:, end+1) = ipVec;
        supThickVec(:, end+1) = thickVec;
        
        % Calculate the ratio and thickness vecs
        supThickVec(:,end+1) = sum([supGcVec(:,end),supIpVec(:,end)],2,'includenan');
        supRatioVec(:,end+1) = supGcVec(:,end)./(supThickVec(:,end) .^ 2);
    end
end
SD = std(supRatioVec .', 'omitnan');
meanSupRatioVecProfile = nanmean(supRatioVec,2);
errorbar(posEccentricity, meanSupRatioVecProfile, SD, '-s','MarkerSize',...
    2, 'MarkerEdgeColor','black','MarkerFaceColor','black')
title('Mean GC:GCIP^2 Along the Superior Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 20])
ylim([0 75])

% plot inferior meridian
infThickVec = [];
infGcVec = [];
infIpVec = [];
infRatioVec = [];

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
        
        infGcVec(:, end+1) = gcVec;
        infIpVec(:, end+1) = ipVec;
        infThickVec(:, end+1) = thickVec;
        
        % Calculate the ratio and thickness vecs
        infThickVec(:,end+1) = sum([infGcVec(:,end),infIpVec(:,end)],2,'includenan');
        infRatioVec(:,end+1) = infGcVec(:,end)./(infThickVec(:,end) .^ 2);
    end
end
SD = std(infRatioVec .', 'omitnan');
meanInfRatioVecProfile = nanmean(infRatioVec,2);
errorbar(posEccentricity, meanInfRatioVecProfile, SD, '-s','MarkerSize',...
    2, 'MarkerEdgeColor','black','MarkerFaceColor','black')
title('Mean GC:GCIP^2 Along the Inferior Meridian')
xlabel('Eccentricity') 
ylabel('GCL:GCIPL^2 Ratio')
xlim([0 20])
ylim([0 75])

%% plot standard deviation across meridans for individual subjects
f4 = figure;

sdProfile = [];

for ii = 1:50
    load(horizThicknessProfile,'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', 'XPos_Degs', 'subIDs');
    posEccentricity = XPos_Degs(1282:2561);

    % temporal meridian
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

        % Calculate the thickness vec and flip vecs
        thickVec = sum([gcVec,ipVec],2,'includenan');
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);

        tempGcVec = gcVec;
        tempIpVec = ipVec;
        tempThickVec = thickVec;
        tempRatioVec = tempGcVec./(tempThickVec .^ 2);
    end

    % nasal meridian
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

        nasGcVec = gcVec;
        nasIpVec = ipVec;
        nasThickVec = thickVec;
        nasRatioVec = nasGcVec(:,end)./(nasThickVec .^ 2);
    end

    % superior meridian
    load(vertThicknessProfile,'GCthicknessValuesAtXPos_um', ...
        'IPthicknessValuesAtXPos_um', 'XPos_Degs');
    posEccentricity = XPos_Degs(1282:2561);

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
        gcVec = flipud(gcVec);
        thickVec = flipud(thickVec);
        ipVec = flipud(ipVec);

        supGcVec = gcVec;
        supIpVec = ipVec;
        supThickVec = thickVec;
        supRatioVec = supGcVec./(supThickVec .^ 2);
    end

    % inferior meridian
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

        infGcVec = gcVec;
        infIpVec = ipVec;
        infThickVec = thickVec;
        infRatioVec = infGcVec./(infThickVec .^ 2);
    end

     subRatioVec = [tempRatioVec, nasRatioVec, supRatioVec, infRatioVec];
     SD = std(subRatioVec .', 'omitnan');
%     meanSubRatioVecProfile = nanmean(subRatioVec, 2);
%     errorbar(posEccentricity, meanSubRatioVecProfile, SD, '-s','MarkerSize',...
%         2, 'MarkerEdgeColor','black','MarkerFaceColor','black')
%     title('Standard Deviation GC:GCIP^2 Across Meridians')
%     xlabel('Eccentricity') 
%     ylabel('GCL:GCIPL^2 Ratio')
%     xlim([0 30])
%     ylim([0 75])
%     hold on
    
    sdProfile(:,end+1) = SD;
end

sdMeanProfile = nanmean(sdProfile, 2);
plot(posEccentricity, sdMeanProfile)
title('Mean Standard Deviation GC:GCIP^2 Across Meridians')
xlabel('Eccentricity') 
ylabel('Standard Deviation')
xlim([0 30])
ylim([0 10])