%% fitModelToAllData
% This routine loads the outputs of processDensityMaps.m and then fits the
% data with the polar cone density surface model.

% The overall result directory
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';

% Identify the aggregate data files
mergedFiles = dir([sourceDir '*_merged.mat']);
subNames = strrep(extractfield(mergedFiles,'name'),'_merged.mat','');

% Define some constants
supportLengthMm = 2001;
maxSupportMm = 5;
supportMmDelta = 0.0025;

% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportMm = 0:supportMmDelta:supportMmDelta*(supportLengthMm-1);


% Load the data from all subjects
dataMatMm = nan(supportLengthMm,supportLengthMm,length(subNames));
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
%    assert(abs(supportMmDelta-data.meta.mmPerPixelFixed)<0.001);

    % Store the polar maps
    dataMatMm(:,:,ss)=data.polarDensityMm(:,:);
    
end

% Filter out any zero or negative values that slipped in
dataMatMm(dataMatMm<=0)=nan;

% Fit the mean mm map
Y = nanmean(dataMatMm,3);
w = sum(~isnan(dataMatMm),3);
p0 = fitDensitySurfaceMm(Y,w,false,false,true,false);

% Fit each subject with the reduced model
pSet = nan(20,length(subNames));
YfitSet = nan(size(dataMatMm));
YResidualSet = nan(size(dataMatMm));
fValSet= nan(1,length(subNames));
RSquaredSet= nan(4,length(subNames));
nonlconSet = nan(1,length(subNames));
polarThetaSet = nan(1,length(subNames));
polarMultiplierSet = nan(1,length(subNames));

fprintf('fitting...');
w1 = ones(size(Y));
for ii = 1:length(subNames)
    Y1 = squeeze(dataMatMm(:,:,ii));
    fprintf([num2str(ii),'...']);
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii), RSquaredSet(:,ii), nonlconSet(ii), polarThetaSet(ii), polarMultiplierSet(ii)] = fitDensitySurfaceMm(Y1,w1,true,true,true,false,p0);
    YResidualSet(:,:,ii) = Y1 - squeeze(YfitSet(:,:,ii));
end
fprintf('done\n');

% Save the individual subject fits
individualFitFile = fullfile(sourceDir,'individualSubjectFitsMm.mat');
save(individualFitFile,'p0','pSet','YfitSet','fValSet','RSquaredSet','polarThetaSet','polarMultiplierSet','dataMatMm','subNames','YResidualSet')


