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
ratiosOD = [];
ratiosOS = [];

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
        
        ratiosOD(end+1) = nanmedian(ratioVecOD);
        ratiosOS(end+1) = nanmedian(ratioVecOS);
    end
end

% Present a plot that demonstrates that there is individual variation in
% the thickness of the GC layer relative to the GC+IP thickness.
plot(ratiosOD,ratiosOS,'xk');
axis square
xlabel('median GC/(GC+IP) ratio OD');
ylabel('median GC/(GC+IP) ratio OD');
refline([1 0]);
corr(ratiosOD',ratiosOS')



end
