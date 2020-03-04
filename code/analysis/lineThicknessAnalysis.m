function comboTable=lineThicknessAnalysis(GCIPthicknessFile, varargin)
% Do some analysis
%
% Description:
%   Foo
%
% Examples:
%{
    dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
    GCIPthicknessFile = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg');
    comboTable=lineThicknessAnalysis(GCIPthicknessFile);
%}

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('GCIPthicknessFile',@ischar);



% Optional analysis params

p.addParameter('showPlots',true,@islogical);
p.addParameter('dataSaveName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'), 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
p.addParameter('mmPerDegFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat'),@ischar);
p.addParameter('figSaveDir','/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/horizontalLine',@ischar);
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

% p.addParameter('dataSaveName',fullfile('C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)', 'AOSO_analysis','OCTExplorerExtendedHorizontalData','LineAnalysisResults.mat'),@ischar);
% p.addParameter('figSaveDir','C:\Users\dontm\Dropbox (Personal)\Research\Publications\Connectome_RetinaAnalysis_2019\figures',@ischar);
% p.addParameter('subjectTableFileName',fullfile('C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)','TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);


%% Parse and check the parameters
p.parse(GCIPthicknessFile, varargin{:});

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% Load the data file
load(GCIPthicknessFile)

% Load the mmPerDegMaps file
load(p.Results.mmPerDegFileName,'mmPerDegPolyFit');

subList = {};
thickVec = [];
gcVec = [];
ipVec = [];
ratioVec = [];
gcMedianOD = [];
gcMedianOS = [];
ipMedianOD = [];
ipMedianOS = [];

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
            gcVec(:,end+1) = nanmean([gcVecOD,gcVecOS],2);
            ipVec(:,end+1) = nanmean([ipVecOD,ipVecOS],2);
        end
        
        % Calculate the ratio and thickness vecs
        thickVec(:,end+1) = sum([gcVec(:,end),ipVec(:,end)],2,'includenan');
        ratioVec(:,end+1) = gcVec(:,end)./thickVec(:,end);
        
        % Save the median value for each eye and layer
        gcMedianOD(end+1) = nanmedian(gcVecOD);
        gcMedianOS(end+1) = nanmedian(gcVecOS);
        ipMedianOD(end+1) = nanmedian(ipVecOD);
        ipMedianOS(end+1) = nanmedian(ipVecOS);
        
    end
end

% Make some vectors of mean thickness and ratio
subCountPerPoint = sum(~isnan(thickVec),2);
meanThickVec = nanmean(thickVec,2);
meanGCVec = nanmean(gcVec,2);
meanRatioVec = nanmean(ratioVec,2);
semThickVec = nanstd(thickVec,1,2)./sqrt(subCountPerPoint);
semRatioVec = nanstd(ratioVec,1,2)./sqrt(subCountPerPoint);
badIdx = subCountPerPoint<(length(subList)/3);
meanGCVec(badIdx)=nan;

% Plot the GC thickness functions
profilePlot(XPos_Degs, gcVec, meanGCVec, 'Eccentricity [deg visual angle]','Thickness [microns]', ...
    ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

% Perform the calculation for the GC layer
totalVolumePerDegSq = zeros(size(gcVec));
gcVolumePerDegSq = zeros(size(gcVec));
AreaPerDegSq = zeros(size(gcVec));

for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    mmSqPerDegSq = mmPerDegPolyFit{idx}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
    axialLengths(ss) = subjectTable.Axial_Length_average(idx);
    totalVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq;
    gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmSqPerDegSq;
end
meanGCVolumePerDegSqProfile = nanmean(gcVolumePerDegSq,2);
meanGCVolumePerDegSqProfile(badIdx) = nan;

% Plot gc tissue volume functions
profilePlot(XPos_Degs, gcVolumePerDegSq, meanGCVolumePerDegSqProfile, 'Eccentricity [deg visual angle]','GC tissue volume [mm^3] / deg^2', ...
    ['GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)


%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmean(gcVec)'),num2cell(nanmean(gcVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcMeanThick','gcVolumePerDegSq'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');

% Plot GC thickness vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.gcMeanThick, 'Axial length [mm]','median GC thickness [microns]', ...
    ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMeanThick))],p.Results.showPlots)

% Plot mean GC tissue volume vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.gcVolumePerDegSq, 'Axial length [mm]','median GC tissue volume [mm^3] / deg^2', ...
    ['Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq))],p.Results.showPlots)

% Slope of axial length vs. mean gcTissueVolume
pp = polyfit(axialLengths,nanmean(gcVolumePerDegSq,1),1);
adjustProportion = polyval(pp,23.58) ./ polyval(pp,axialLengths) ;

for ss = 1:50
    gcVolumePerDegSqAdjust(:,ss) = gcVolumePerDegSq(:,ss) .* adjustProportion(ss);
end

meanGCVolumePerDegSqProfileAdjust = nanmean(gcVolumePerDegSqAdjust,2);
meanGCVolumePerDegSqProfileAdjust(badIdx) = nan;

% Plot adjusted gc tissue volume functions
profilePlot(XPos_Degs, gcVolumePerDegSqAdjust, meanGCVolumePerDegSqProfileAdjust, 'Eccentricity [deg visual angle]','GC tissue volume [mm^3] / deg^2', ...
    ['Adjusted GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)


% Plot adjusted mean GC tissue volume vs axial length.
regressionPlot(axialLengths, nanmean(gcVolumePerDegSqAdjust,1), 'Axial length [mm]','adjusted mean GC tissue volume [mm^3] / deg^2', ...
    ['Axial length vs. adjusted gc tissue volume, r=',num2str(corr(axialLengths',nanmean(gcVolumePerDegSqAdjust,1)'))],p.Results.showPlots)


%% PCA

% Set to nan any x positions for which more than 10 subjects are missing
% data
gcVolumePerDegSqAdjust(badIdx,:) = nan;
nanX = sum(isnan(gcVolumePerDegSqAdjust'))>47;
gcVolumeCleaned = gcVolumePerDegSqAdjust(~nanX,:);

% Conduct a PCA using the alternating least squares (ALS) algorithm to
% account for the missing values
[coeff,score,~,~,explained,mu] = pca(gcVolumeCleaned,'Centered',false,'algorithm','als');

% Expand the score vectors back to the full x range
scoreExpanded = nan(size(gcVolumePerDegSqAdjust));
scoreExpanded(~nanX,:) = score;

% Find the three segment domains
tt = find(diff(nanX));
domains = {[1 tt(1)],[tt(2)+1 tt(3)],[tt(4)+1 length(nanX)]};

% Perform piece-wise spline smoothing of the scores to remove the noisy
% effects of data imputation
smoothVal = 0.1; % 0-1, lower is smoother.
scoreExpandedSmoothed = scoreExpanded;
for cc = 1:size(scoreExpanded,2)
    for dd = 1:length(domains)
        rd = domains{dd};
        x = XPos_Degs(rd(1):rd(2));
        y = scoreExpanded(rd(1):rd(2),cc);
        pp = csaps(x',y,smoothVal);
        scoreExpandedSmoothed(rd(1):rd(2),cc) = ppval(pp,x');
    end
end

% Show the effect of smoothing on the components
nDimsToUse = 6;
figure
set(gcf,'color','w');
for ii = 1:nDimsToUse
subplot(3,2,ii);
plot(scoreExpanded(:,ii),'.','Color',[0.85 0.85 0.85]);
hold on
plot(scoreExpandedSmoothed(:,ii),'-r','LineWidth',1);
axis off
end

% Plot the reconstructions
figure
set(gcf,'color','w');
for ii = 1:49
subplot(7,7,ii);
plot(gcVolumePerDegSqAdjust(:,ii),'.','Color',[0.85 0.85 0.85]);
hold on
profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff(ii,1:nDimsToUse)';
plot(profileFit,'-r','LineWidth',1);
axis off
end

% Plot the PCs
profilePlot(XPos_Degs, scoreExpandedSmoothed(:,1).*(coeff(:,1)'), scoreExpandedSmoothed(:,1).*(mean(coeff(:,1))), 'Eccentricity [deg visual angle]','PC1', ...
    ['PC1 for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)
xlim([-25 25]);
hold on
plot(XPos_Degs,meanGCVolumePerDegSqProfileAdjust,'-g');


end

