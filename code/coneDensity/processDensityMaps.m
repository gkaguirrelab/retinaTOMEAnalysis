% processDensityMaps
%
% This routine loads the output of Rob Cooper's cone density estimation
% routines, and then combines the "merged" data (composed of the confocal
% and split images) and the "foveal" single image. Along the way there are
% some cleaning steps applied, including filtering outlier density values
% and removing an artifactual decline in cone density close to the fovea in
% the confocal data.
%
% The images are all aligned to the estimated foveal center, and
% transformed to polar coordinates.
%
% The code is written to support a separate processing of a "merged" and
% "fovea" dataset. Rob ultimately added the foveal images into the
% aggregate image, so this aspect of the code is vestigal.
%

clear

% Hard coded values
downSample = 0.05; % Downsamples the image to 1/20th of the original rez
newDim = 18000; % Dimensions of the density maps
pixelsPerDegreeFixed = 647; % All of the results files have this rez.
paraFovealExtent = 2; % Start of search for "ridge" to filter (in degrees)

% Create a map that will be used to filter "extreme" values
dParams = [1477 -0.3396 7846 -1.3049 629];
density = @(x) dParams(1).*exp(dParams(2).*x)+dParams(3).*exp(dParams(4).*x)+dParams(5);
[R,C] = ndgrid(1:newDim, 1:newDim);
R = R - newDim/2;
C = C - newDim/2;
r = sqrt(R.^2+C.^2) ./ pixelsPerDegreeFixed;
r(r<1.25)=nan;

% The overal result directory
dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));

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

% We are going to loop through the processing once for the aggregate
% images. There is vestigal code here that would have handled a separate
% set of foveal images, but these are now integrated into the "aggregate"
% images.
for dd = 1:1

    switch dd
        case 1
            resultFiles=dir('Aggregation_Analysis/*_merged_Fouriest_Result.mat');
            tagName = '_merged';
            maxThresh = density(r).*2.5;
        case 2
            resultFiles=dir('*/confocal/Results_fovea_map/*_Fovea.mat');
            tagName = '_fovea';
            maxThresh = density(r).*4;
    end

    % Loop over the result files
    for rr = 1:length(resultFiles)

        % Identify the next file
        fileName = fullfile(resultFiles(rr).folder,resultFiles(rr).name);

        % Extract the subject name
        switch dd
            case 1
                tmp = strsplit(fileName,filesep);
                subName = [tmp{end}(4:8) tmp{end}(18:20)];
            case 2
                tmp = strsplit(fileName,filesep);
                subName = tmp{end-3};
        end

        % Report that we are about to process this subject
        fprintf([resultFiles(rr).name '\n']);

        % Load the file
        load(fileName);

        % Make sure that the image resolution is as expected
        if dd==1
            assert(1/scaling==pixelsPerDegreeFixed)
        end

        % The density map information is stored in different places in the
        % merged and "fovea" datasets.
        switch dd
            case 1
                density_map = density_map_comb;
            case 2
                density_map = densim(:,:,1);
        end

        % Determine if this a left or right eye
        if contains(resultFiles(rr).name,'_OS_')
            laterality = 'OS';
            polarToMeridian = {'Inferior','Temporal','Superior','Nasal'};
        elseif contains(resultFiles(rr).name,'_OD_')
            laterality = 'OD';
            polarToMeridian = {'Inferior','Nasal','Superior','Temporal'};
        else
            error(['Unable to determine laterality for ' resultFiles(rr).name]);
        end

        % Code to allow an over-ride of the default fovea coords; not used.
        foveaOverride = false;
        switch dd
            case 1
                % Structure of a special case foveal coordinate
                if strcmp(subName,'123456_OD')
                    fovea_coords = [4444, 5555];
                    foveaOverride = true;
                end
                foveaCoordStore.(['s_' subName])=fovea_coords;
            case 2
                fovea_coords=foveaCoordStore.(['s_' subName]);
        end

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

        % In the fovea image, filter to find just the highest connected
        % components
        if dd == 2
            imDensity = findOrderedRegions(imDensity);
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

        % Store the data
        data = [];
        data.meta.subName = subName;
        data.meta.tagName = tagName;
        data.meta.filename = resultFiles(rr).name;
        data.meta.folder = resultFiles(rr).folder;
        data.meta.laterality = laterality;
        data.meta.downSample = downSample;
        data.meta.foveaCoords = fovea_coords;
        data.meta.foveaOverride = foveaOverride;
        data.meta.supportDegDelta = supportDeg(1);
        data.imDensity = imDensity;
        data.polarDensity = polarDensity;

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

end

% Restore the warning state
warning(warnState);

