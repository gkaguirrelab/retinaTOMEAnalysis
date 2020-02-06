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

subList = {};
thickVec = [];
gcVec = [];
ipVec = [];
ratioVec = [];
gcMedianOD = [];
gcMedianOS = [];
ipMedianOD = [];
ipMedianOS = [];

% Obtain the GC+IP thickness and ratio functions for each subject. While we
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
badIdx = subCountPerPoint<(length(subList)/2);
meanThickVec(badIdx)=nan;
meanGCVec(badIdx)=nan;
meanRatioVec(badIdx)=nan;
semThickVec(badIdx)=nan;
semRatioVec(badIdx)=nan;


close all
% Plot the GC thickness functions
profilePlot(XPos_Degs, gcVec, meanGCVec, 'Eccentricity [deg visual angle]','Thickness [microns]', ...
    ['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

% Plot the GC+IP thickness functions
profilePlot(XPos_Degs, thickVec, meanThickVec, 'Eccentricity [deg visual angle]','Thickness [microns]', ...
    ['GC+IP thickness profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

% Plot the ratio functions
profilePlot(XPos_Degs, ratioVec, meanRatioVec, 'Eccentricity [deg visual angle]','Ratio GC/[GC+IP]', ...
    ['GC ratio profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

% median GC vs. IP thickness
aa = nanmean([gcMedianOD; gcMedianOS]);
bb = nanmean([ipMedianOD; ipMedianOS]);
regressionPlot(aa, bb, 'median GC thickness','median IP thickness', ...
    ['GC vs IP median thickness across subjects, r=',num2str(corr(aa',bb'))],p.Results.showPlots)


% Nasal vs. temporal ratio differences
if p.Results.showPlots
    figure
    subVec = 1:floor(length(XPos_Degs)/2);
    
    % Plot the temporal arm
    plot(XPos_Degs(subVec),meanRatioVec(subVec),'-k','LineWidth',2);
    hold on
    plot(XPos_Degs(subVec),meanRatioVec(subVec)+semRatioVec(subVec),'-k','LineWidth',1);
    plot(XPos_Degs(subVec),meanRatioVec(subVec)-semRatioVec(subVec),'-k','LineWidth',1);
    
    % Now mirror the vectors and plot the nasal arm
    meanRatioVecFlip = flipud(meanRatioVec);
    semRatioVecFlip = flipud(semRatioVec);
    plot(XPos_Degs(subVec),meanRatioVecFlip(subVec),'-r','LineWidth',2);
    plot(XPos_Degs(subVec),meanRatioVecFlip(subVec)+semRatioVecFlip(subVec),'-r','LineWidth',1);
    plot(XPos_Degs(subVec),meanRatioVecFlip(subVec)-semRatioVecFlip(subVec),'-r','LineWidth',1);
    
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Ratio GC/[GC+IP]');
    title(['Mean GC ratio profiles [+-SEM] for nasal (red) and temporal (black)'])
end

% Convert thickness to tissue volume. Start with an anonymous function that
% provides mmPerDeg at the ellipsoidal pole of the vitreous chamber
mmPerDeg = @(axialLength) (0.0165.*axialLength)-0.1070;

% Perform the calculation for the GC layer
totalVolumePerDegSq = zeros(size(gcVec));
gcVolumePerDegSq = zeros(size(gcVec));
AreaPerDegSq = zeros(size(gcVec));

for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    axialLength = subjectTable.Axial_Length_average(idx);
    totalVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmPerDeg(axialLength).^2;
    gcVolumePerDegSq(:,ss) = (gcVec(:,ss)).*mmPerDeg(axialLength).^2;
    AreaPerDegSq(:,ss) = mmPerDeg(axialLength).^2;
end
meanGCVolumePerDegSq = nanmean(gcVolumePerDegSq,2);
meanGCVolumePerDegSq(badIdx) = nan;

% Perform the calculation for the IP layer
ipVolumePerDegSq = zeros(size(ipVec));
for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    axialLength = subjectTable.Axial_Length_average(idx);
    ipVolumePerDegSq(:,ss) = ipVec(:,ss).*mmPerDeg(axialLength).^2;
end
meanIPVolumePerDegSq = nanmean(ipVolumePerDegSq,2);
meanIPVolumePerDegSq(badIdx) = nan;


% Plot gc tissue volume functions
profilePlot(XPos_Degs, gcVolumePerDegSq, meanGCVolumePerDegSq, 'Eccentricity [deg visual angle]','GC tissue volume [mm^3] / deg^2', ...
    ['GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))],p.Results.showPlots)

%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmedian(gcVec)'),num2cell(nanmedian(ipVec)'),num2cell(nanmedian(thickVec)'),num2cell(nanmedian(ratioVec)'),num2cell(nanmedian(gcVolumePerDegSq)'),num2cell(nanmedian(ipVolumePerDegSq)')],...
    'VariableNames',{'AOSO_ID','gcMedianThick','ipMedianThick','gcipMedianThick','medianRatio','gcVolumePerDegSq','ipVolumePerDegSq'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');

% Plot GC thickness vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.gcMedianThick, 'Axial length [mm]','median GC thickness [microns]', ...
    ['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMedianThick))],p.Results.showPlots)

% Plot IP thickness vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.ipMedianThick, 'Axial length [mm]','median IP thickness [microns]', ...
    ['Axial length vs. median IP thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.ipMedianThick))],p.Results.showPlots)

% Plot ratio vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.medianRatio, 'Axial length [mm]','median GC/(GC+IP) ratio', ...
    ['Axial length vs. median ratio, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.medianRatio))],p.Results.showPlots)

% Plot median GC tissue volume vs axial length.
regressionPlot(comboTable.Axial_Length_average, comboTable.gcVolumePerDegSq, 'Axial length [mm]','median GC tissue volume [mm^3] / deg^2', ...
    ['Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq))],p.Results.showPlots)

% Plot median GC thickness vs median GC tissue volume.
regressionPlot(comboTable.gcMedianThick, comboTable.gcVolumePerDegSq, 'median GC thickness [micron]','median GC tissue volume [mm^3] / deg^2', ...
    ['gc thickness vs. gc tissue volume, r=',num2str(corr(comboTable.gcMedianThick,comboTable.gcVolumePerDegSq))],p.Results.showPlots)


% Plot median GC+IP tissue volume vs axial length.
sumVec = comboTable.gcVolumePerDegSq+comboTable.ipVolumePerDegSq;
regressionPlot(comboTable.Axial_Length_average, sumVec, 'Axial length [mm]','median GC+IP tissue volume [mm^3] / deg^2', ...
    ['Axial length vs. gc+ip tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,sumVec))],p.Results.showPlots)

%p.Results.showPlots=1;
% Plot r values of GC tissue volume vs axial length across Ecc
if p.Results.showPlots
    figure
    corrGCVolumePerDegSq=zeros(1,size(gcVolumePerDegSq,1));
    for i = 1:size(gcVolumePerDegSq,1)
        corrGCVolumePerDegSq(i)= corr(comboTable.Axial_Length_average,gcVolumePerDegSq(i,:)', 'rows','complete');
    end
    plot(XPos_Degs,corrGCVolumePerDegSq,'-k','LineWidth',4);
    xlabel('Eccentricity (deg)');
    ylabel('Correlation Coefficient, r ');
end


% Plot linear slope of GC tissue volume vs axial length across Ecc
if p.Results.showPlots
    figure
    [YData, YData_errTop,YData_errBot] = slopeAtEcc(gcVolumePerDegSq,comboTable.Axial_Length_average,0,0);
    
    plot(XPos_Degs,YData,'-k','LineWidth',2);
    xlabel('Eccentricity (deg)');
    ylabel('Regression Slope');
    %    ylabel('Regression Intercept');
    %    ylabel('Regression Slope As Percent of Total Volume');
    
    figure
    plot(XPos_Degs,YData,'-k','LineWidth',2);
    hold on
    plot(XPos_Degs,squeeze(YData_errTop),'-r','LineWidth',2);
    plot(XPos_Degs,squeeze(YData_errBot),'-r','LineWidth',2);
    hold off
    
    xlabel('Eccentricity (deg)');
    ylabel('Regression Slope');
    %    ylabel('Regression Intercept');
    %    ylabel('Regression Slope As Percent of Total Volume');
    
end

% Plot linear slope of area expansion vs axial length across Ecc
if p.Results.showPlots
    figure
    [YData2, YData_errTop,YData_errBot] = slopeAtEcc(AreaPerDegSq,comboTable.Axial_Length_average,0,0);
    
    plot(XPos_Degs,YData2,'-k','LineWidth',2);
    xlabel('Eccentricity (deg)');
    %    ylabel('Regression Slope');
    %    ylabel('Regression Intercept');
    ylabel('Regression Slope');
    
    figure
    plot(XPos_Degs,YData2,'-k','LineWidth',2);
    hold on
    plot(XPos_Degs,squeeze(YData_errTop),'-r','LineWidth',2);
    plot(XPos_Degs,squeeze(YData_errBot),'-r','LineWidth',2);
    hold off
    xlabel('Eccentricity (deg)');
    ylabel('Regression Slope');
    
    figure
    plot(XPos_Degs,YData./YData2,'-k','LineWidth',2);
    
end

%plot slope of gc volume vs area
if p.Results.showPlots
    figure
    [YData2, YData_errTop,YData_errBot] = slopeAtEcc2(gcVolumePerDegSq,AreaPerDegSq,0,0);
    
    plot(XPos_Degs,YData2,'-k','LineWidth',2);
    xlabel('Eccentricity (deg)');
    %    ylabel('Regression Slope');
    %    ylabel('Regression Intercept');
    ylabel('Regression Slope');
    
    figure
    plot(XPos_Degs,YData2,'-k','LineWidth',2);
    hold on
    plot(XPos_Degs,squeeze(YData_errTop),'-r','LineWidth',2);
    plot(XPos_Degs,squeeze(YData_errBot),'-r','LineWidth',2);
    hold off
    xlabel('Eccentricity (deg)');
    ylabel('Regression Slope');
    
end


%PCA Anlaysis
%work in progress
% Need to deal with Nans here
%k=0:.05:30;
k=0:.01:1;%step size
w=nanmean(nanmean(gcVec,2))';
slope = zeros(1,length(k));
for i = 1:length(k)
slope(i) = gcPCAcorrCostFunction(k(i),totalVolumePerDegSq,gcVec,AreaPerDegSq,comboTable.Axial_Length_average,XPos_Degs);
end

%plot the cose function
figure
plot(k,abs(slope))
axis square
xlabel('Stromal Component Thickness (fraction of PC1 of GC profile) ');
ylabel('Slope of Median GC Volume vs Axial Length');

%find the optimal point of the cost function
[~,minInd]=min(abs(slope));

[slope_optimal, adjustedgcVolumePerDegSq_optimal] = gcPCAcorrCostFunction(k(minInd),totalVolumePerDegSq,gcVec,AreaPerDegSq,comboTable.Axial_Length_average,XPos_Degs);


%plot regression at current point
adjustedgcVolumePerDegSq_optimal_median=nanmedian(adjustedgcVolumePerDegSq_optimal,2);
regressionPlot(comboTable.Axial_Length_average, adjustedgcVolumePerDegSq_optimal_median, 'Axial length [mm]','Adjusted tissue volume [mm^3] / deg^2', ...
    ['Axial length vs. adjusted gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,adjustedgcVolumePerDegSq_optimal_median))],p.Results.showPlots)

%plot volume profile of each subject after adjusment
profilePlot(XPos_Degs, adjustedgcVolumePerDegSq_optimal, nanmean(adjustedgcVolumePerDegSq_optimal,1), 'Eccentricity [deg visual angle]','Adjusted tissue volume [mm^3] / deg^2', ...
    ['GC Volume profiles for each subject (and mean) after adjustment, n=',num2str(length(subList))],1)

%example of one profile before and after
figure
plot(XPos_Degs, adjustedgcVolumePerDegSq_optimal(10,:));
hold on
plot(XPos_Degs, gcVolumePerDegSq(:,10));


save(p.Results.dataSaveName,'XPos_Degs','meanGCVec');
end
