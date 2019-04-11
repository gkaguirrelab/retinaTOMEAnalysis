function analyzeMps(thicknessMapDir,volumeMapDir, saveDir, varargin)
% Do some analysis
%
% Description:
%   Foo
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('thicknessMapDir',@ischar);
p.addRequired('volumeMapDir',@ischar);
p.addRequired('saveDir',@ischar);

% Optional analysis params
p.addParameter('degreesFOV',30,@isscalar);
p.addParameter('showPlots',false,@islogical);
p.addParameter('layerSetLabels',{'RGCIPL','RNFL','OPL','TotalRetina'},@iscell);

%% Parse and check the parameters
p.parse(thicknessMapDir, volumeMapDir, saveDir, varargin{:});

subIDs = dir(fullfile(volumeMapDir,'1*'));

%We're going to look at some measurements, this describes
measurements = zeros(length(subIDs),9);%this will be out ouput matrix
%This describes each column in the measurements
header = {'subject ID', 'RGCIPL mean thickness', 'RGCIPL mean volume', ...
    'RNFL mean thickness', 'RNFL mean volume', ...
    'OPL mean thickness', 'OPL mean volume', ...
    'Total Retina mean thickness', 'Total Retina mean volume'};


for layer = 1:length(p.Results.layerSetLabels) %L controls which layer we're looking
    overlap = [];
    
    %first we find the overlaping region across all subjects
    for ss = 1:length(subIDs)
        LoadthicknessMap=load(fullfile(thicknessMapDir,subIDs(ss).name,[subIDs(ss).name '_averageMaps.mat']));
        
        thicknessMap = LoadthicknessMap.averageMaps.(p.Results.layerSetLabels{layer});
        loadname = fullfile(volumeMapDir, subIDs(ss).name, [subIDs(ss).name '_' p.Results.layerSetLabels{layer} '_volumeMap.mat']);
        load(loadname,'volumeMap_mmCubedDegSquared');
        
        if(isempty(overlap))
            overlap = ones(size(volumeMap_mmCubedDegSquared));
        end
        
        overlap= overlap & ~isnan(volumeMap_mmCubedDegSquared) &  ~isnan(thicknessMap);
    end
    
    % writeout the overlap
    
    savename = fullfile(saveDir, [p.Results.layerSetLabels{layer} ' _overlapMap.mat']);
    save(savename,'overlap');
    
    
    %now that we've got the overlap, we go back through and calculate the
    %mean for each subject across the overlap
    
    for ss = 1:length(subIDs)
        LoadthicknessMap=load(fullfile(thicknessMapDir,subIDs(ss).name,[subIDs(ss).name '_averageMaps.mat']));
        
        thicknessMap = LoadthicknessMap.averageMaps.(p.Results.layerSetLabels{layer});
        loadname = fullfile(volumeMapDir, subIDs(ss).name, [subIDs(ss).name '_' p.Results.layerSetLabels{layer} '_volumeMap.mat']);
        load(loadname,'volumeMap_mmCubedDegSquared');
        
        measurements(ss,1) = str2double(subIDs(ss).name);
        measurements(ss,2*layer) = mean(thicknessMap(overlap));
        measurements(ss,2*layer+1) =  mean(volumeMap_mmCubedDegSquared(overlap));
    end
    
end


filename = fullfile(saveDir,'meanThicknessAndVolumes.xlsx');
xlswrite(filename,header,1,'A1');
xlswrite(filename,measurements,1,'A2');