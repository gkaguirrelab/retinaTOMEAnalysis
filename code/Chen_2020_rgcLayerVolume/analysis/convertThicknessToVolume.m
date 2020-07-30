function [mmSqPerDegSq,gcVolumePerDegSq,meanGCVolumePerDegSqProfile,volumeTable] = convertThicknessToVolume(p,gcVec,badIdx,subList,XPos_Degs,subjectTable,orientation)
%% Convert from mm thickness to tissue volume

% Define some variables
gcVolumePerDegSq = zeros(size(gcVec));

% Load the mmPerDegMaps file
load(p.Results.mmPerDegFileName,'mmPerDegPolyFit');


% Loop over subjects
for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    switch orientation
        case 'horiz'
            mmSqPerDegSq(:,ss) = mmPerDegPolyFit{idx}([-XPos_Degs;zeros(size(XPos_Degs))]').^2;
        case 'vert'
            mmSqPerDegSq(:,ss) = mmPerDegPolyFit{idx}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
    end
    gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);
end

% Get the mean tissue volume profile, and nan out the "bad" indices
meanGCVolumePerDegSqProfile = nanmean(gcVolumePerDegSq,2);
meanGCVolumePerDegSqProfile(badIdx) = nan;

volumeTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcVolumePerDegSq'});