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



% Plot the GC thickness functions
if p.Results.showPlots
    figure
    plot(XPos_Degs,gcVec,'-r');
    hold on
    plot(XPos_Degs,meanGCVec,'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Thickness [microns]');
    title(['GC thickness profiles for each subject (and mean), n=',num2str(length(subList))])
end

% Plot the GC+IP thickness functions
if p.Results.showPlots
    figure
    plot(XPos_Degs,thickVec,'-r');
    hold on
    plot(XPos_Degs,meanThickVec,'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Thickness [microns]');
    title(['GC+IP thickness profiles for each subject (and mean), n=',num2str(length(subList))])
end

% Plot the ratio functions
if p.Results.showPlots
    figure
    plot(XPos_Degs,ratioVec,'-r');
    hold on
    plot(XPos_Degs,meanRatioVec,'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Ratio GC/[GC+IP]');
    title(['GC ratio profiles for each subject (and mean), n=',num2str(length(subList))])
end

% median GC vs. IP thickness
if p.Results.showPlots
    figure
    aa = nanmean([gcMedianOD; gcMedianOS]);
    bb = nanmean([ipMedianOD; ipMedianOS]);    
    plot(aa,bb,'xr');
    hold on
    c = polyfit(aa,bb,1);
    plot(aa,polyval(c,aa),'--b')

    axis square
    xlabel('median GC thickness');
    ylabel('median IP thickness');
    title(['GC vs IP median thickness across subjects, r=',num2str(corr(aa',bb'))])
end

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
gcVolumePerDegSq = zeros(size(gcVec));
for ss = 1:length(subList)
    idx = find(subjectTable.AOSO_ID==str2num(subList{ss}));
    axialLength = subjectTable.Axial_Length_average(idx);
    gcVolumePerDegSq(:,ss) = gcVec(:,ss).*mmPerDeg(axialLength).^2;
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
if p.Results.showPlots
    figure
    plot(XPos_Degs,gcVolumePerDegSq,'-r');
    hold on
    plot(XPos_Degs,meanGCVolumePerDegSq,'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('GC tissue volume [mm^3] / deg^2');
    title(['GC tissue volume profiles for each subject (and mean), n=',num2str(length(subList))])
end


%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmedian(gcVec)'),num2cell(nanmedian(ipVec)'),num2cell(nanmedian(thickVec)'),num2cell(nanmedian(ratioVec)'),num2cell(nanmedian(gcVolumePerDegSq)'),num2cell(nanmedian(ipVolumePerDegSq)')],'VariableNames',{'AOSO_ID','gcMedianThick','ipMedianThick','gcipMedianThick','medianRatio','gcVolumePerDegSq','ipVolumePerDegSq'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');

% Plot GC thickness vs axial length.
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.gcMedianThick,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.gcMedianThick,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC thickness [microns]');
    title(['Axial length vs. median GC thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcMedianThick))])
end

% Plot IP thickness vs axial length.
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.ipMedianThick,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.ipMedianThick,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median IP thickness [microns]');
    title(['Axial length vs. median IP thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.ipMedianThick))])
end


% Plot ratio vs axial length.
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.medianRatio,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.medianRatio,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC/(GC+IP) ratio');
    title(['Axial length vs. median ratio, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.medianRatio))])
end

% Plot median GC tissue volume vs axial length.
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC tissue volume [mm^3] / deg^2');
    title(['Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq))])
end

% Need to deal with Nans here
[coeff,score,~,~,explained,mu] = pca(gcVolumePerDegSq,'Centered',true);
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC tissue volume [mm^3] / deg^2');
    title(['Axial length vs. gc tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcVolumePerDegSq))])
end


% Plot median GC+IP tissue volume vs axial length.
if p.Results.showPlots
    figure
    sumVec = comboTable.gcVolumePerDegSq+comboTable.ipVolumePerDegSq;
    plot(comboTable.Axial_Length_average,sumVec,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,sumVec,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC+IP tissue volume [mm^3] / deg^2');
    title(['Axial length vs. gc+ip tissue volume, r=',num2str(corr(comboTable.Axial_Length_average,sumVec))])
end

save(p.Results.dataSaveName,'XPos_Degs','meanGCVec');

end
