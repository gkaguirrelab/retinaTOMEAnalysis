% processDensityMaps
%
% This routine loads the output of Rob Cooper's cone density estimation
% routines, which produce a cone density map from a combination of
% confocal, split, and a single confocal foveal image. Along the way there
% are some cleaning steps applied, including filtering outlier density
% values and removing an artifactual decline in cone density close to the
% fovea.
%
% The images are all aligned to the estimated foveal center, and
% transformed to polar coordinates.
%

clear

% Hard coded values
downSample = 0.05; % Downsamples the image to 1/20th of the original rez
newDim = 18000; % Dimensions of the density maps
pixelsPerDegreeFixed = 647; % All of the results files have this rez.
paraFovealExtent = 2; % Start of search for "ridge" to filter (in degrees)
mmPerPixelFixed = 0.01; % Map resolution in retinal mm coordinate space

% Create a map that will be used to filter "extreme" values
dParams = [1477 -0.3396 7846 -1.3049 629];
density = @(x) dParams(1).*exp(dParams(2).*x)+dParams(3).*exp(dParams(4).*x)+dParams(5);
[R,C] = ndgrid(1:newDim, 1:newDim);
R = R - newDim/2;
C = C - newDim/2;
r = sqrt(R.^2+C.^2) ./ pixelsPerDegreeFixed;
r(r<1.25)=nan;
maxThresh = density(r).*2.5;

% The overal result directory
dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));

% Load the subject data table
subjectTableFileName = fullfile(dropboxBaseDir,'TOME_subject','TOME-AOSO_SubjectInfo.xlsx');
opts = detectImportOptions(subjectTableFileName);
subjectTable = readtable(subjectTableFileName, opts);

% Change to our working directory
cd(fullfile(dropboxBaseDir,'Connectome_AOmontages_images'))

% Load the foveaCoordsStore if it exists
fileName = fullfile('densityAnalysis','foveaCoordStore.mat');
if isfile(fileName)
    load(fileName,'foveaCoordStore');
end

% Turn off a warning
warnState = warning();
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off','signal:findpeaks:largeMinPeakHeight');

resultFiles=dir('Aggregation_Analysis/*_merged_Fouriest_Result.mat');
tagName = '_merged';

% Loop over the result files
for rr = 1:length(resultFiles)

    % Identify the next file
    fileName = fullfile(resultFiles(rr).folder,resultFiles(rr).name);

    % Extract the subject name
    tmp = strsplit(fileName,filesep);
    tmp = tmp{end};
    tmpBits = strsplit(tmp,'_');
    subName = tmpBits{2}; subDate = tmpBits{3}; subEye = tmpBits{4};

    % Report that we are about to process this subject
    fprintf([subName '_' subEye '\n']);

    % Make sure that this subject has a row in the table
    subRow = find(subjectTable.AOSO_ID==str2double(subName));

    if length(subRow) ~= 1
        warning('This subject is missing an entry in the subject table; skipping')
        continue
    end

    % Load the file
    load(fileName);

    % Make sure that the image resolution is as expected
    assert(1/scaling==pixelsPerDegreeFixed)

    % Grab the density map
    density_map = density_map_comb;

    % Determine if this a left or right eye
    switch subEye
        case 'OS'
            laterality = 'OS';
            polarToMeridian = {'Inferior','Temporal','Superior','Nasal'};
        case 'OD'
            laterality = 'OD';
            polarToMeridian = {'Inferior','Nasal','Superior','Temporal'};
        otherwise
            error(['Unable to determine laterality for ' resultFiles(rr).name]);
    end

    % Code to allow an over-ride of the default fovea coords
    foveaOverride = false;
    if strcmp(subName,'11098') && strcmp(subEye,'OS')
        fovea_coords = [7025, 8013];
        foveaOverride = true;
    end
    foveaCoordStore.(['s_' subName])=fovea_coords;

    % Pad the density map so that it is square around the fovea_coords
    imsize = size(density_map);
    offset = round(newDim/2-fovea_coords);

    % Set up the density map
    imDensity = zeros(newDim,newDim);
    imDensity(1:imsize(1),1:imsize(2))=density_map;
    imDensity(isnan(imDensity))=0;
    imDensity(isinf(imDensity))=0;
    imDensity=imtranslate(imDensity,offset);
    imDensity(imDensity==0)=nan;

    % Show the initial density map
    figHandle = figure('Name',subName,'Position',[100, 100, 800, 200]);
    subplot(1,3,1)
    imagesc(imDensity)
    hold on
    plot([1 newDim],round([newDim newDim]/2),'-r');
    plot(round([newDim newDim]/2),[1 newDim],'-r');
    axis square
    axis off
    title('centered')

    % Filter the map for extreme values
    imDensity(imDensity>maxThresh)=nan;

    % Down-sample the map
    imDensity = imresize(imDensity,downSample);

    % Flip left eye density maps so they are pseudo-right eye
    panelTitle = 'resample, filter';
    if strcmp(laterality,'OS')
        imDensity = fliplr(imDensity);
        panelTitle = 'flip, resample, filter';
    end

    % Check if imDensity is empty, in which case throw an error
    if sum(~isnan(imDensity(:)))==0
        error('No data points survived!')
    end

    % Show the filtered map
    subplot(1,3,2)
    imagesc(imDensity)
    hold on
    plot([1 newDim],round([newDim newDim]/2),'-r');
    plot(round([newDim newDim]/2),[1 newDim],'-r');
    axis square
    axis off
    title(panelTitle)

    % Create polar maps
    polarDensity = convertImageMapToPolarMap(imDensity);
    polarDim = newDim*downSample*2-1;

    % Show the polar density map pre filtering
    subplot(1,3,3)
    imagesc(polarDensity)
    axis square
    axis off
    title('polar w/ ridge filter')

    % Calculate the support in degrees for the polar image
    supportDeg = (1:polarDim)./(pixelsPerDegreeFixed*downSample*4);

    % Remove decreasing components close to the fovea
    supportDegIdx = find(supportDeg>paraFovealExtent,1);
    threshIdx = ones(1,size(polarDensity,1));

    % Loop over the polar angles to find the ridge
    for nn=1:size(polarDensity,1)

        % Density for this polar angle
        myVec = polarDensity(nn,1:supportDegIdx);

        [~,idx] = findpeaks(myVec,'MinPeakHeight',4e3,'MinPeakProminence',500,'MinPeakWidth',10,'NPeaks',3);

        % Couldn't find a peak. Try this method instead
        if isempty(idx)
            stillSearching = true;
            idx = 1;
            while stillSearching
                if myVec(idx) > max(myVec(idx+1:end))
                    stillSearching = false;
                else
                    idx=idx+1;
                    if idx==length(myVec)
                        stillSearching = false;
                    end
                end
            end
        end

        % Store the ridge index
        threshIdx(nn)=idx(end)-1;

    end % loop over polar angle

    % Draw the filtering ridge
    hold on
    plot(threshIdx,1:size(polarDensity,1),'-r');
    drawnow

    % Apply the filtering
    for nn=1:size(polarDensity,1)
        polarDensity(nn,1:threshIdx(nn))=nan;
    end

    % Convert the cleaned polar image back to Cartesian
    imDensity = convertPolarMapToImageMap(polarDensity);

    % Check if imDensity is empty, in which case throw an error
    if sum(~isnan(imDensity(:)))==0
        error('No data points survived!')
    end

    % Create a model eye, and a warp field to bring the data into mm
    % retinal coordinates
    axialLength = subjectTable.Axial_Length_average(subRow);
    sphericalAmetropia = subjectTable.Spherical_Error_average(subRow);
    eye = modelEyeParameters('axialLength',axialLength,'sphericalAmetropia',sphericalAmetropia,'accommodation',0);

    % Define some landmarks needed for visual field calculation
    rayOriginDistance = min([2000, 1000 / eye.meta.accommodation]);
    principalPoint = calcPrincipalPoint(eye, rayOriginDistance);

    % Define the parameters of a displacement map. The warp field is quite
    % smooth, so we define the map at 0.5 mm increments, out to 5 mm away
    % from the fovea
    imRes = 0.5;
    imExtent = 5;
    dim = (imExtent*2)/imRes+1;
    D = nan(dim,dim,2);
    xVals = ((1:dim)-ceil(dim/2))*imRes;
    yVals = ((1:dim)-ceil(dim/2))*imRes;
    for xx = 1:dim
        for yy = 1:dim
            theta = wrapTo360(atan2d(yVals(yy),xVals(xx)));
            eccen = sqrt(xVals(xx)^2+yVals(yy)^2);

            X = calcRetina2DPolToCart(eye,theta,eccen);
            [~,fieldAngularPosition] = calcNodalRayToRetina(eye,X,rayOriginDistance,principalPoint);
            fieldAngularPosition = fieldAngularPosition - eye.landmarks.fovea.degField;

            D(xx,yy,:) = fieldAngularPosition;
        end
    end

    % Slight numerical errors cause the foveal coordinate to not have a
    % visual field position of exactly zero. Fudge that here
    D(:,:,1) = D(:,:,1) - D(ceil(dim/2),ceil(dim/2),1);
    D(:,:,2) = D(:,:,2) - D(ceil(dim/2),ceil(dim/2),2);

    % Interpolate the warp field to the full, desired resolution of 10
    % microns.
    dimHi = (imExtent*2)/mmPerPixelFixed+1;
    xValsHi = ((1:dimHi)-ceil(dimHi/2))*mmPerPixelFixed;
    yValsHi = ((1:dimHi)-ceil(dimHi/2))*mmPerPixelFixed;
    DHi(:,:,1) = interp2(xVals,yVals',squeeze(D(:,:,1)),xValsHi,yValsHi');
    DHi(:,:,2) = interp2(xVals,yVals',squeeze(D(:,:,2)),xValsHi,yValsHi');

    % Adjust the values to map pixels from one space to the other
    DWarp = (DHi ./ (supportDeg(1)*4)) + ceil(size(imDensity,1)/2);
    [H,V]=meshgrid(1:dimHi,1:dimHi);
    DWarp(:,:,1) = squeeze(DWarp(:,:,1)) - H;
    DWarp(:,:,2) = squeeze(DWarp(:,:,2)) - V;

    % Get the density map in mm coordinates
    imDensityMm = imwarp(imDensity,DWarp);

    % Convert to polar coordinates
    polarDensityMm = convertImageMapToPolarMap(imDensityMm);

    % Store the meta
    data = [];
    data.meta.subName = subName;
    data.meta.subDate = subDate;
    data.meta.subEye = subEye;
    data.meta.tagName = tagName;
    data.meta.filename = resultFiles(rr).name;
    data.meta.folder = resultFiles(rr).folder;
    data.meta.laterality = laterality;
    data.meta.downSample = downSample;
    data.meta.foveaCoords = fovea_coords;
    data.meta.foveaOverride = foveaOverride;
    data.meta.supportDegDelta = supportDeg(1);
    data.meta.mmPerPixelFixed = mmPerPixelFixed;
    data.meta.axialLength = axialLength;
    data.meta.sphericalAmetropia = sphericalAmetropia;

    % Store the data
    data.imDensity = imDensity;
    data.polarDensity = polarDensity;
    data.imDensityMm = imDensityMm;
    data.polarDensityMm = polarDensityMm;
    data.DWarp = DWarp;
    data.eye = eye;

    % Save the data file
    fileName = fullfile('densityAnalysis',[subName tagName '.mat']);
    save(fileName,'data','-v7.3');

    % Save the foveaCoordsStore
    fileName = fullfile('densityAnalysis','foveaCoordStore.mat');
    save(fileName,'foveaCoordStore');

    % Save the diagnostic image
    fileName = fullfile('densityAnalysis',[subName tagName '.png']);
    saveas(figHandle,fileName)

    % Close the image
    close(figHandle);

    % Clear the data file
    clear data

end


% Restore the warning state
warning(warnState);

