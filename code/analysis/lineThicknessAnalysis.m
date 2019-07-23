function lineThicknessAnalysis(GCIPthicknessFile, varargin)
% Do some analysis
%
% Description:
%   Foo
%
% Examples:
%{
    dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
    GCIPthicknessFile = fullfile(dropboxBaseDir, 'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg_7_18_2019.mat');
    lineThicknessAnalysis(GCIPthicknessFile)
%}

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('dataDir',@ischar);

% Optional analysis params
p.addParameter('showPlots',true,@islogical);
p.addParameter('figSaveDir','/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/horizontalLine',@ischar);
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

%% Parse and check the parameters
p.parse(GCIPthicknessFile, varargin{:});

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);
axialLength = subjectTable.Axial_Length_average;

% Load the data file
GCIPthicknessFile =    '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/OCTExplorerExtendedHorizontalData/GCIP_thicknessesByDeg_7_18_2019.mat';
load(GCIPthicknessFile)

subList = {};
thickVec = [];
ratioVec = [];
ratioOD = [];
ratioOS = [];
thickOD = [];
thickOS = [];

% Obtain the GC+IP thickness and ratio functions for each subject. While we
% are at it, confirm that there is a substantial correlation across
% subjects between the left and right eye in the median of the ratio
% functions.
GCthicknessValuesAtXPos(GCthicknessValuesAtXPos==0)=nan;
IPthicknessValuesAtXPos(IPthicknessValuesAtXPos==0)=nan;

for ii = 1:50
    if ~all(isnan(squeeze(GCthicknessValuesAtXPos(ii,1,:)))) && ...
            ~all(isnan(squeeze(GCthicknessValuesAtXPos(ii,2,:))))
        gcVecOD = squeeze(GCthicknessValuesAtXPos(ii,1,:));
        gcVecOS = flipud(squeeze(GCthicknessValuesAtXPos(ii,2,:)));
        gcVec = mean([gcVecOD,gcVecOS],2,'includenan');
        ipVecOD = squeeze(IPthicknessValuesAtXPos(ii,1,:));
        ipVecOS = flipud(squeeze(IPthicknessValuesAtXPos(ii,2,:)));
        ipVec = mean([ipVecOD,ipVecOS],2,'includenan');
        
        thickVecOD = sum([gcVecOD,ipVecOD],2,'includenan');
        thickVecOS = sum([gcVecOS,ipVecOS],2,'includenan');
        thickVec(:,end+1) = sum([gcVec,ipVec],2,'includenan');
        ratioVecOD = gcVecOD./thickVecOD;
        ratioVecOS = gcVecOS./thickVecOS;
        ratioVec(:,end+1) = gcVec./thickVec(:,end);
        subList(end+1) = {subIDs(ii,:)};
        
        % Save some values from the left and right eye for diagnostic plots
        thickOD(end+1) = nanmedian(thickVecOD);
        thickOS(end+1) = nanmedian(thickVecOS);
        ratioOD(end+1) = nanmedian(ratioVecOD);
        ratioOS(end+1) = nanmedian(ratioVecOS);
    end
end

% Present a plot that demonstrates that there is individual variation in
% the thickness of the GC+IP layer.
if p.Results.showPlots
    figure
    plot(thickOD,thickOS,'xr');
    axis square
    xlabel('median GC+IP thickness [pixels] OD');
    ylabel('median GC+IP thickness [pixels] OS');
    title(['Individual variation in median GC+IP thickness, r=',num2str(corr(thickOD',thickOS'))])
    refline([1 0]);
end

% Present a plot that demonstrates that there is individual variation in
% the thickness of the GC layer relative to the GC+IP thickness.
if p.Results.showPlots
    figure
    plot(ratioOD,ratioOS,'xr');
    axis square
    xlabel('median GC/(GC+IP) ratio OD');
    ylabel('median GC/(GC+IP) ratio OS');
    title(['Individual variation in median GC/[GC+IP] thickness, r=',num2str(corr(ratioOD',ratioOS'))])
    refline([1 0]);
end

% Plot the thickness functions
if p.Results.showPlots
    figure
    plot(XPos_Degs,thickVec,'-r');
    hold on
    plot(XPos_Degs,nanmean(thickVec,2),'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Thickness [microns]');
    title(['GC+IP thickness profiles for each subject (and mean), n=',num2str(length(subList))])
end

% Plot the ratio functions
if p.Results.showPlots
    figure
    plot(XPos_Degs,ratioVec,'-r');
    hold on
    plot(XPos_Degs,nanmean(ratioVec,2),'-k','LineWidth',4);
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Ratio GC/[GC+IP]');
    title(['GC ratio profiles for each subject (and mean), n=',num2str(length(subList))])
end

% median ratio vs. thickness
if p.Results.showPlots
    figure
    plot(nanmedian(thickVec),nanmedian(ratioVec),'xr');
    hold on
    c = polyfit(nanmedian(thickVec),nanmedian(ratioVec),1);
    plot(nanmedian(thickVec),polyval(c,nanmedian(thickVec)),'--b')

    axis square
    xlabel('median GC+IP thickness');
    ylabel('median GC/(GC+IP) ratio');
    title(['Individual variation in median GC/[GC+IP] thickness, r=',num2str(corr(nanmedian(thickVec)',nanmedian(ratioVec)'))])
end

% Nasal vs. temporal ratio differences
if p.Results.showPlots
    figure
    ratioMean = nanmean(ratioVec,2);
    ratioSEM = nanstd(ratioVec,1,2)./sqrt(length(subList));
    subVec = 1:floor(length(XPos_Degs)/2);

    % Plot the temporal arm
    plot(XPos_Degs(subVec),ratioMean(subVec),'-k','LineWidth',2);
    hold on
    plot(XPos_Degs(subVec),ratioMean(subVec)+ratioSEM(subVec),'-k','LineWidth',1);
    plot(XPos_Degs(subVec),ratioMean(subVec)-ratioSEM(subVec),'-k','LineWidth',1);

    % Now mirror the vectors and plot the nasal arm
    ratioMean = flipud(ratioMean);
    ratioSEM = flipud(ratioSEM);
    plot(XPos_Degs(subVec),ratioMean(subVec),'-r','LineWidth',2);
    plot(XPos_Degs(subVec),ratioMean(subVec)+ratioSEM(subVec),'-r','LineWidth',1);
    plot(XPos_Degs(subVec),ratioMean(subVec)-ratioSEM(subVec),'-r','LineWidth',1);
    
    xlabel('Eccentricity [deg visual angle]');
    ylabel('Ratio GC/[GC+IP]');
    title(['Mean GC ratio profiles [+-SEM] for nasal (red) and temporal (black)'])
end


%% Relate GC+IP thickness and ratio to axial length

% Create a table of median thickness and axial length
dataTable = cell2table([num2cell(str2double(subList)'),num2cell(nanmedian(thickVec)'),num2cell(nanmedian(ratioVec)')],'VariableNames',{'AOSO_ID','gcipMedianThick','medianRatio'});

% Join the data table with the subject biometry and demographics table,
% using the AOSO_ID as the key variable
comboTable = join(dataTable,subjectTable,'Keys','AOSO_ID');

% Plot GC+IP thickness vs axial length.
if p.Results.showPlots
    figure
    plot(comboTable.Axial_Length_average,comboTable.gcipMedianThick,'xr');
    hold on
    c = polyfit(comboTable.Axial_Length_average,comboTable.gcipMedianThick,1);
    plot(comboTable.Axial_Length_average,polyval(c,comboTable.Axial_Length_average),'--b')
    axis square
    xlabel('Axial length [mm]');
    ylabel('median GC+IP thickness [microns]');
    title(['Axial length vs. median GC+IP thickness, r=',num2str(corr(comboTable.Axial_Length_average,comboTable.gcipMedianThick))])
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


end
