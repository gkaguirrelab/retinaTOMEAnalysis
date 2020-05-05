%% Convert from mm thickness to tissue volume

% Define some variables
totalVolumePerDegSq = zeros(size(gcVec));
gcVolumePerDegSq = zeros(size(gcVec));

% Loop over subjects
for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    mmSqPerDegSq(:,ss) = mmPerDegPolyFit{idx}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
    axialLengths(ss) = subjectTable.Axial_Length_average(idx);
    totalVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);
    gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);
end

% Get the mean tissue volume profile, and nan out the "bad" indices
meanGCVolumePerDegSqProfile = nanmean(gcVolumePerDegSq,2);
meanGCVolumePerDegSqProfile(badIdx) = nan;
