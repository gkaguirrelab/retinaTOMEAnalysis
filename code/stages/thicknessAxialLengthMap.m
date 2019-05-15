function thicknessAxialLengthMap(dataDir, varargin)
% Do some analysis
%
% Description:
%   Foo
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('dataDir',@ischar);

% Optional analysis params
p.addParameter('layerSetLabels',{'RGCIPL','RNFL'},@iscell);
p.addParameter('showPlots',true,@islogical);
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

%% Parse and check the parameters
p.parse(dataDir, varargin{:});


% Obtain a list of subjects
rawSubjectList = dir(fullfile(dataDir,'*/*.mat'));
nSubs = length(rawSubjectList);


% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);
axialLength = subjectTable.Axial_Length_average;

% Plot
if p.Results.showPlots
    figure
end
    
% Loop over layer sets
for ii = 1:length(p.Results.layerSetLabels)
    
    % Loop over subjects and load the maps
    for ss = 1:length(rawSubjectList)
        fileName = fullfile(rawSubjectList(ss).folder,rawSubjectList(ss).name);
        load(fileName,'averageMaps');
        thisMap = averageMaps.(p.Results.layerSetLabels{ii});
        if ss==1
            imageSize = size(thisMap);
            avgMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
        end
        avgMapBySubject(ss,:,:)=thisMap;
    end
    
    % Loop over positions in the map and calculate the correlation between
    % axial length and layer thickness across subjects
    corrMap = nan(imageSize(1),imageSize(2));
    for xx=1:imageSize(1)
        for yy=1:imageSize(2)
            corrMap(xx,yy) = corr2(squeeze(avgMapBySubject(:,xx,yy)),axialLength);
        end
    end
    
    % Plot
    if p.Results.showPlots
        subplot(2,1,ii)
            mesh(corrMap);
            caxis([-1 1]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([-1 1]);
            title(p.Results.layerSetLabels{ii})
            axis square
    end
    
end
