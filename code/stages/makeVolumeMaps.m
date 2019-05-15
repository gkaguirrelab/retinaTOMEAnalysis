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
p.addParameter('layerSetLabels',{'RGCIPL'},@iscell);

%% Parse and check the parameters
p.parse(thicknessMapDir, mmPerDegMapDir, saveDir, varargin{:});


octRadialDegreesVisualExtent = p.Results.degreesFOV/2;

%find all subjects
subIDs = dir(fullfile(p.Results.thicknessMapDir,'1*'));

% Loop over subjects
for ii = 1:length(subIDs)

    % Load average thickness maps
    load(fullfile(p.Results.thicknessMapDir,subIDs(ii).name,[subIDs(ii).name '_averageMaps.mat']),'averageMaps');

    % Load the mmPerDeg map
    load(fullfile(p.Results.mmPerDegMapDir,[subIDs(ii).name '_mmPerDegMap.mat']),'mmPerDeg');
    
    % We will fit a thin plate spline to the mmPerDeg map to allow
    % interpolation. First define the support in visual degrees for the
    % mmPerDeg map
    supportDegLowRes = linspace(-octRadialDegreesVisualExtent,octRadialDegreesVisualExtent,size(mmPerDeg,1));

    % Silence a warning that occurs regarding nans in the maps
    warningState = warning;
    warning('off','curvefit:prepareFittingData:removingNaNAndInf');

    % Prepare the surface anf perform the fit
    [xo,yo,zo]=prepareSurfaceData(supportDegLowRes,supportDegLowRes,mmPerDeg);
    mmPerDegFit = fit([xo, yo],zo,'thinplateinterp');

    % Restore the warning state
    warning(warningState);

    % Create an interpolated mmPerDeg
    dimX = size(averageMaps.(p.Results.layerSetLabels{1}),1);
    dimY = size(averageMaps.(p.Results.layerSetLabels{1}),2);
    supportDegHiRes = linspace(-octRadialDegreesVisualExtent,octRadialDegreesVisualExtent,dimX);
    [xi,yi]=meshgrid(supportDegHiRes,supportDegHiRes);
    mmPerDegHiRes = mmPerDegFit(xi,yi);
    
    % Create a directory for the output for this subject
    outDir = fullfile(p.Results.saveDir, subIDs(ii).name);
    if ~exist(outDir)
        mkdir(outDir);
    end
    
    % Loop over the layer sets that we wish to process
    for jj = 1:length(p.Results.layerSetLabels)
        
        % Obtain the thickness map for this layer set
        thicknessMicronMap = averageMaps.(p.Results.layerSetLabels{jj});
                
        % Re-orient the mmPerDeg map to the OD clinical view
        mmPerDegHiResRot = fliplr(mmPerDegHiRes');
        
        % Calculate the degrees per pixel
        degreesPerPixel = p.Results.degreesFOV / dimX;
        
        % Convert thickness map to volume map with thickness(mm)*degree^2
        % The value at each point in the resulting map is the point-wise
        % estimate of the volume (in mm) of retinal tissue per deg^2 of
        % visual field.
        volumeMap_mmCubedDegSquared = (thicknessMicronMap./1000).*(mmPerDegHiResRot.^2);
        
        % Save the volume map
        savename = fullfile(p.Results.saveDir, subIDs(ii).name, [subIDs(ii).name '_' p.Results.layerSetLabels{jj} '_volumeMap.mat']);
        save(savename,'volumeMap_mmCubedDegSquared');
        
    end
        
end
