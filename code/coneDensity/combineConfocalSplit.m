
% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')

% Load the confocal and split data
dataFileName = 'confocalDensityProfileData.mat';
load(dataFileName)
confocalData = data;

dataFileName = 'splitDensityProfileData.mat';
load(dataFileName)
splitData = data;

clear data


%% Aggregate the profiles
% Find the longest support deg
[confocalSupportLength,idx] = max(cellfun(@(x) length(x.profile.supportDeg),confocalData));
confocalSupportDeg = 0:confocalData{1}.meta.supportDegDelta:confocalData{1}.meta.supportDegDelta*(confocalSupportLength-1);

[splitSupportLength,idx] = max(cellfun(@(x) length(x.profile.supportDeg),splitData));
splitSupportDeg = 0:splitData{1}.meta.supportDegDelta:splitData{1}.meta.supportDegDelta*(splitSupportLength-1);

polarToMeridian = {'Inferior','Nasal','Superior','Temporal'};
polarToVal = [270,0,90,180];

figure
confocalStartIdx = find(confocalSupportDeg > 0.5,1);
splitStartIdx = find(splitSupportDeg > 0.5,1);

% Loop over the arms
for mm = 1:length(polarToMeridian)
    confocalDataMatrix = nan(length(confocalData),confocalSupportLength);
    for rr = 1:length(confocalData)
        tmp = confocalData{rr}.profile.(polarToMeridian{mm});
        confocalDataMatrix(rr,1:length(tmp))=tmp;
    end
    
    splitDataMatrix = nan(length(confocalData),splitSupportLength);
    for rr = 1:length(splitData)
        tmp = splitData{rr}.profile.(polarToMeridian{mm});
        splitDataMatrix(rr,1:length(tmp))=tmp;
    end
    
    % Add the lines to the plot
    subplot(2,2,mm);
    tmp = nanmedian(confocalDataMatrix);
    plot(confocalSupportDeg(confocalStartIdx:end),tmp(confocalStartIdx:end),'-');
    
    hold on
    tmp = nanmedian(splitDataMatrix);
    plot(splitSupportDeg(splitStartIdx:end),tmp(splitStartIdx:end),'-');
    
%    fitHandle = getSplineFitToConeDensitySqDegVisual(polarToVal(mm));
%    plot(splitSupportDeg,fitHandle(splitSupportDeg),':k');
    
    title(polarToMeridian{mm});
    ylim([0 7000]);
    xlim([0 15]);
    ylabel('mean density [cones/deg^2]');
    xlabel('distance from fovea [deg]')
    
end
