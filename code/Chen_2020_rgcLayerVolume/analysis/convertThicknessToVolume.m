function [mmSqPerDegSq,gcVolumePerDegSq,meanGCVolumePerDegSqProfile,volumeTable] = convertThicknessToVolume(p,gcVec,badIdx,subList,XPos_Degs,subjectTable,orientation)
%% Convert from mm thickness to tissue volume

% Define some variables
gcVolumePerDegSq = zeros(size(gcVec));

% Load the eye models file
load(p.Results.eyeModelsFileName,'eyeModels');

% Create the down-sampled XPos_Degs that we will measure
XPos_Support = linspace(1,length(XPos_Degs),100);
XPos_DegSub = interp1(1:length(XPos_Degs),XPos_Degs,XPos_Support);

% Loop over subjects
for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    eye = eyeModels{idx};
    fieldAngularPosition = eye.landmarks.fovea.degField;
    rayOriginDistance = 1500;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    mmPerDeg = nan(size(XPos_DegSub));
    for ii = 1:length(XPos_DegSub)
        switch orientation
            case 'horiz'
                deltaAngles = [1 0]; % Measure horizontally separated points
                thisPosition = fieldAngularPosition + [ -XPos_DegSub(ii) 0 ];
            case 'vert'
                deltaAngles = [1 0]; % Measure vertically separated points
                thisPosition = fieldAngularPosition + [ 0 -XPos_DegSub(ii) ];
        end
        mmPerDeg(ii) = calcMmRetinaPerDeg(eye,thisPosition,deltaAngles,rayOriginDistance,angleReferenceCoord);
    end
    mmSqPerDegSq(:,ss) = interp1(XPos_DegSub,mmPerDeg,XPos_Degs).^2;
    gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq(:,ss);
end

% Get the mean tissue volume profile, and nan out the "bad" indices
meanGCVolumePerDegSqProfile = nanmean(gcVolumePerDegSq,2);
meanGCVolumePerDegSqProfile(badIdx) = nan;

volumeTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcVolumePerDegSq'});