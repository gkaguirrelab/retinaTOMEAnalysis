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



%% Parse and check the parameters
p.parse(GCIPthicknessFile, varargin{:});

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% Load the data file
load(GCIPthicknessFile,'XPos_Degs', ...
    'GCthicknessValuesAtXPos_um', ...
    'IPthicknessValuesAtXPos_um', ...
    'subIDs');

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

% Define the "bad" indices across x position as those that are missing
% measurements from more than a 2/3rds of the subjects.
badIdx = subCountPerPoint<(length(subList)/3);
meanGCVec(badIdx)=nan;

% Plot the GC thickness functions
profilePlot(XPos_Degs, gcVec, meanGCVec, 'Eccentricity [deg visual angle]','Thickness [microns]', ...
    ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)


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

% Plot gc tissue volume profile
profilePlot(XPos_Degs, gcVolumePerDegSq, meanGCVolumePerDegSqProfile, 'Eccentricity [deg visual angle]','GC tissue volume [mm^3] / deg^2', ...
    ['GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

% Plot the mmSqPerDegSq functions
profilePlot(XPos_Degs, mmSqPerDegSq, nanmean(mmSqPerDegSq,2), 'Eccentricity [deg visual angle]','mm^2 / deg^2', ...
    ['mm square retina per degree square visuall angle (and mean), n=',num2str(length(subList))],p.Results.showPlots)


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



%% Conduct PCA upon the tissue volume data

% Before we do the PCA, nan out the "bad" indices from the data. We won't
% attempt to reconstruct these points.
gcVolumePerDegSq(badIdx,:) = nan;

% We limit the PCA analysis to those x positions for which we have
% measurements for greater than 90% of the subjects.
nanX = sum(isnan(gcVolumePerDegSq'))>45;
gcVolumeCleaned = gcVolumePerDegSq(~nanX,:);

% Conduct a PCA using the alternating least squares (ALS) algorithm to
% account for the few missing values
[coeff,score,~,~,explained] = pca(gcVolumeCleaned,'Centered',false,'algorithm','als');

% Expand the score vectors back to the full x range
scoreExpanded = nan(size(gcVolumePerDegSq));
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

% Show the effect of smoothing on the PCA scores
nDimsToUse = 6;
figure
set(gcf,'color','w');
shapeExplained = explained ./ sum(explained(2:end));
shapeExplained(1) = nan;
for ii = 1:nDimsToUse
    subplot(3,2,ii);
    plot(scoreExpanded(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    plot(scoreExpandedSmoothed(:,ii),'-r','LineWidth',1);
    str = sprintf('PC%d, shape var explained: %2.2f',ii,shapeExplained(ii));
    title(str);
    axis off
end
suptitle('Original and smoothed PCA scores')

% Plot the reconstructions
figure
set(gcf,'color','w');
for ii = 1:49
    subplot(7,7,ii);
    plot(gcVolumePerDegSq(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff(ii,1:nDimsToUse)';
    plot(profileFit,'-r','LineWidth',1);
    axis off
end
suptitle('Original and fitted gc tissue volume profiles by subject')



%% Identify the scores that correlate with axial length
corrCoeffAxialLength = nan(1,nDimsToUse);
fprintf('Correlation axial length with PCA coefficients:\n');
for ii = 1:nDimsToUse
    corrCoeffAxialLength(ii) = corr(axialLengths',coeff(:,ii));
    txt = sprintf('\tPC%d: %2.2f \n',ii,corrCoeffAxialLength(ii));
    fprintf(txt);
end

% Regress the axial length upon the coefficients
X = [coeff(:,1:nDimsToUse)'; ones(1,50)]';
b = X\axialLengths';

% Show regression fit
figure
plot(axialLengths,X*b,'*r');
axis square
refline(1,0);

%%%%%%%%%%%%%%%
%%
%% NEED TO ORTHOGONALIZE X WITH RESPECT TO AXIAL LENGTH HERE
%%
%% AND ADD BACK IN THE EFFECT OF THE EMMETROPIC EYE (23.58 mm)
%%
%% also synthesize
adjustedCoeff = zeros(size(coeff));
ALmax= max(axialLengths);
ALmin= min(axialLengths);
ALRange = linspace(ALmin,ALmax,50);
ALRange(16) = 23.58;%We're forcing this one to be the mean emmetropic value (was 23.56 before).
synCoeff = zeros(50,nDimsToUse);
for p=1:nDimsToUse
    
    currCoeff = coeff(:,p);
    %find regression line between PC coeff and  axial length
    pp = polyfit(axialLengths,currCoeff',1);
    
    %find PC coeff value for emmetrope
    %Adjust each score by fraction: coeff_emm / coeff_regress
    AL_emmetrope_mm=23.58;
    adjustProportion_allSub = polyval(pp,AL_emmetrope_mm) - polyval(pp,axialLengths);
    adjustedCoeff(:,p) = currCoeff + adjustProportion_allSub';
    
    %Save synthesized representation of AL
    synCoeff(:,p) = polyval(pp,ALRange);
end

%%%%%%%%%%%%%%%
% Plot the reconstructions with the adjustment
figure
set(gcf,'color','w');
[ALsorted, ALsortedIndx] = sort(axialLengths);
counter=0;
for ii = ALsortedIndx
    counter = counter+1;
    subplot(8,7,counter);
    plot(gcVolumePerDegSq(:,ii),'.','Color',[0.85 0.85 0.85]);
    hold on
    profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*adjustedCoeff(ii,1:nDimsToUse)';
    plot(profileFit,'-r','LineWidth',1);
    ylim([-1 8])
    axis off
end
suptitle('Original and AL influence adjusted gc tissue volume profiles by subject')

%%%%%%%%%%%%%%%
% Plot the synthesized reconstructions by axial length
figure
set(gcf,'color','w');
for ii = 1:50
    subplot(8,7,ii);
    profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(ii,1:nDimsToUse)';
    plot(profileFit,'-r','LineWidth',1);
    ylim([-1 8])
    axis off
end
suptitle('Synthesized gc tissue volume profiles by AL')

figure

profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(1,1:nDimsToUse)';
plot(profileFit,'-r','LineWidth',1);
hold on
profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(16,1:nDimsToUse)';
plot(profileFit,'-g','LineWidth',1);
profileFit = scoreExpandedSmoothed(:,1:nDimsToUse)*synCoeff(50,1:nDimsToUse)';
plot(profileFit,'-b','LineWidth',1);
ylim([-1 8])
title('Synthetic from Axial Lengths: Red=21.79, Green=23.58, Blue=27.57')

%%plot each PC by axial length
figure

%use center profile if we want to see it relative to an emmetropic eye
%(23.58mm)
centerProfile = scoreExpandedSmoothed(:,1)*synCoeff(16,1)';

for d = 2:6
subplot(3,2,d)
profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(1,d)';
plot(profileFit,'-r','LineWidth',1);
hold on
profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(16,d)';
plot(profileFit,'-g','LineWidth',1);
profileFit = centerProfile + scoreExpandedSmoothed(:,d)*synCoeff(50,d)';
plot(profileFit,'-b','LineWidth',1);
end
suptitle('PC components by Axial Lengths: Red=21.79, Green=23.58, Blue=27.57')


%
% % Plot the PCs
% meanPC1 = scoreExpandedSmoothed(:,1).*(mean(coeff(:,1)));
% profilePlot(XPos_Degs, scoreExpandedSmoothed(:,1).*(coeff(:,1)'), meanPC1, 'Eccentricity [deg visual angle]','PC1', ...
%     ['PC1 for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)
% xlim([-25 25]);
% hold on
% plot(XPos_Degs,meanGCVolumePerDegSqProfileAdjust,'-g');

% %% Adjust gcTissue volume to emmetropic eye
%
% % Slope of axial length vs. mean gcTissueVolume
% pp = polyfit(axialLengths,nanmean(gcVolumePerDegSq,1),1);
% adjustProportion = polyval(pp,23.58) ./ polyval(pp,axialLengths) ;
%
% % Loop through the subjects and apply the adjust proportion to the profile
% for ss = 1:50
%     gcVolumePerDegSqAdjust(:,ss) = gcVolumePerDegSq(:,ss) .* adjustProportion(ss);
% end
%
% % Obtain the mean, adjusted profile, and nan out the "bad" indices.
% meanGCVolumePerDegSqProfileAdjust = nanmean(gcVolumePerDegSqAdjust,2);
% meanGCVolumePerDegSqProfileAdjust(badIdx) = nan;
%
% % Plot adjusted gc tissue volume functions
% profilePlot(XPos_Degs, gcVolumePerDegSqAdjust, meanGCVolumePerDegSqProfileAdjust, 'Eccentricity [deg visual angle]','GC tissue volume [mm^3] / deg^2', ...
%     ['Adjusted GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)
%
%
% % Plot adjusted mean GC tissue volume vs axial length.
% regressionPlot(axialLengths, nanmean(gcVolumePerDegSqAdjust,1), 'Axial length [mm]','adjusted mean GC tissue volume [mm^3] / deg^2', ...
%     ['Axial length vs. adjusted gc tissue volume, r=',num2str(corr(axialLengths',nanmean(gcVolumePerDegSqAdjust,1)'))],p.Results.showPlots)



%% Save some results
%save(p.Results.dataSaveName,'XPos_Degs','meanPC1','gcVolumePerDegSqAdjust','coeff','scoreExpandedSmoothed');

end

