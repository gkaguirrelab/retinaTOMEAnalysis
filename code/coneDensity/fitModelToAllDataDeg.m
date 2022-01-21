%% fitModelToAllDataDeg
% This routine loads the outputs of processDensityMaps.m and then fits the
% data with the polar cone density surface model.

% The overall result directory
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';

% Identify the aggregate data files
mergedFiles = dir([sourceDir '*_merged.mat']);
subNames = strrep(extractfield(mergedFiles,'name'),'_merged.mat','');

% Define some constants
supportLengthDeg = 1799;
maxSupportDeg = 15;
supportDegDelta = 0.00773;

% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLengthDeg-1);

%% Loop through subjects and create the composite polar density image

dataMatDeg = nan(supportLengthDeg,supportLengthDeg,length(subNames));
missingMerged = false(length(subNames));

for ss = 1:length(subNames)
    
    % Load the aggregate data file
    merFile = fullfile(sourceDir,[subNames{ss} '_merged.mat']);
    if isfile(merFile)
        load(merFile,'data');
    else
        str = fprintf(['No merged data for ' subNames{ss} ]);
        error(str);
    end        
    
    % Make sure that the support is as expected
    assert(abs(supportDegDelta-data.meta.supportDegDelta)<0.001);

    % Store the polar maps
    dataMatDeg(:,:,ss)=data.polarDensity(:,:);
    
end

% Filter out any zero or negative values that slipped in
dataMatDeg(dataMatDeg<=0)=nan;

% Fit the mean degree map
Y = nanmean(dataMatDeg,3);
w = sum(~isnan(dataMatDeg),3);
[p0, Yfit, fVal] = fitDensitySurfaceDeg(Y,w,false,false,true,false);

%% Fit each subject with the reduced model
pSet = nan(20,length(subNames));
YfitSet = nan(size(dataMatDeg));
YResidualSet = nan(size(dataMatDeg));
fValSet= nan(1,length(subNames));
RSquaredSet= nan(4,length(subNames));
nonlconSet = nan(1,length(subNames));
polarThetaSet = nan(1,length(subNames));
polarMultiplierSet = nan(1,length(subNames));

fprintf('fitting...');
w1 = ones(size(Y));
for ii = 1:length(subNames)
    Y1 = squeeze(dataMatDeg(:,:,ii));
    fprintf([num2str(ii),'...']);
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii), RSquaredSet(:,ii), nonlconSet(ii), polarThetaSet(ii), polarMultiplierSet(ii)] = fitDensitySurfaceDeg(Y1,w1,true,true,true,false,p0);
    YResidualSet(:,:,ii) = Y1 - squeeze(YfitSet(:,:,ii));
end
fprintf('done\n');

% Save the individual subject fits
individualFitFile = fullfile(sourceDir,'individualSubjectFitsDeg.mat');
save(individualFitFile,'p0','pSet','YfitSet','fValSet','RSquaredSet','polarThetaSet','polarMultiplierSet','dataMatDeg','subNames','YResidualSet')


