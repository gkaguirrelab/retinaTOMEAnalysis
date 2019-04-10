function makeVolumeMaps(thicknessMapDir,mmPerDegMapDir,saveDir, varargin)
% Creates and saves maps of tissue volume per square degree visual field
%
% Description:
%   Foo
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('thicknessMapDir',@ischar);
p.addRequired('mmPerDegMapDir',@ischar);
p.addRequired('saveDir',@ischar);

% Optional analysis params
p.addParameter('degreesFOV',30,@isscalar);
p.addParameter('showPlots',false,@islogical);


% Optional analysis params
p.addParameter('layerSetLabels',{'RGCIPL','RNFL','OPL','TotalRetina'},@iscell);

%% Parse and check the parameters
p.parse(thicknessMapDir, mmPerDegMapDir, saveDir, varargin{:});


%find all subjects
subIDs = dir(fullfile(p.Results.thicknessMapDir,'1*'));

% Loop over subjects
for ii = 1:length(subIDs)
    
    %load mmPerDeg Maps
    LoadmmPerDeg=load(fullfile(p.Results.mmPerDegMapDir,[subIDs(ii).name '_mmPerDegMap.mat']));
    mmPerDeg=LoadmmPerDeg.mmPerDeg;
    
    %load average thickness maps
    LoadthicknessMap=load(fullfile(p.Results.thicknessMapDir,subIDs(ii).name,[subIDs(ii).name '_averageMaps.mat']));
    
    % Create a directory for the output for this subject
    outDir = fullfile(p.Results.saveDir, subIDs(ii).name);
    if ~exist(outDir)
        mkdir(outDir);
    end
    
    % Loop over the layer sets that we wish to process
    for L = 1:length(p.Results.layerSetLabels)
        thicknessMicronMap = LoadthicknessMap.averageMaps.(p.Results.layerSetLabels{L});
        savename = fullfile(p.Results.saveDir, subIDs(ii).name, [subIDs(ii).name '_' p.Results.layerSetLabels{L} '_volumeMap.mat']);
        
        XN = size(thicknessMicronMap,1);
        YN = size(thicknessMicronMap,2);
        
        %fill in nans
        mmPerDegInterp = fillmissing(mmPerDeg,'linear');
        
        %re-orient image to OD clinical view
        mmPerDegInterp_rot = fliplr(mmPerDegInterp');
        
        %resize to same size as slo
        mmPerDegMapInterp_rot_resize = imresize(mmPerDegInterp_rot,[XN YN]);
        
        % Calculate the degrees per pixel
        degreesPerPixel = p.Results.degreesFOV / XN;
        
        %convert thickness map to volume map with thickness(mm)*degree^2
        volumeMap_mmCubedDegSquared = (thicknessMicronMap./100).*((mmPerDegMapInterp_rot_resize*degreesPerPixel).^2);
        
        
        save(savename,'volumeMap_mmCubedDegSquared');
        
    end

    if p.Results.showPlots
    figure(1)
    imshow(mmPerDegMapInterp_rot_resize)
    caxis([min(mmPerDegMapInterp_rot_resize(:)) max(mmPerDegMapInterp_rot_resize(:))])
    colorbar
    figure(2)
    imshow(volumeMap_mmDegSquared)
    caxis([min(volumeMap_mmDegSquared(:)) max(volumeMap_mmDegSquared(:))])
    colorbar
    figure(3)
    imshow(volumeMap_mmCubed)
    caxis([min(volumeMap_mmCubed(:)) max(volumeMap_mmCubed(:))])
    colorbar
    end
    
end
