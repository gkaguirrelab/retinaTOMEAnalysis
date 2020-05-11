function [gcVec,meanGCVecProfile,badIdx,subList,XPos_Degs,subjectTable,thicknessTable] = loadAndMergeData(p,GCIPthicknessFile)

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% Load the data file
load(GCIPthicknessFile,'XPos_Degs', ...
    'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', ...
    'subIDs');


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
        
        % Get the data for each layer and eye
        gcVecOD = squeeze(GCthicknessValuesAtXPos_um(ii,1,:));
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos_um(ii,2,:)));
        ipVecOD = squeeze(IPthicknessValuesAtXPos_um(ii,1,:));
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos_um(ii,2,:)));
        
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

% Report the correlation of mean thickness between eyes across subjects
str = sprintf('The correlation of mean GCIPL thickness between eyes across subjects is R = %2.2f \n',corr(gcipMeanOD',gcipMeanOS','Rows','complete'));
fprintf(str);

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanThickVecProfile = nanmean(thickVec,2);
meanGCVecProfile = nanmean(gcVec,2);
meanRatioVecProfile = nanmean(ratioVec,2);
semThickVecProfile = nanstd(thickVec,1,2)./sqrt(subCountPerPoint);
semRatioVecProfile = nanstd(ratioVec,1,2)./sqrt(subCountPerPoint);

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanGCVecProfile(badIdx)=nan;

% Create a table of median thickness and axial length
thicknessTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVec)')],...
    'VariableNames',{'AOSO_ID','gcMeanThick'});

