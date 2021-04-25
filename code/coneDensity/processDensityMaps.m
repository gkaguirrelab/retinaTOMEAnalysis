clear

% Hard coded values
downSample = 0.05;
foveaDilate = 100;
newDim = 18000; % Dimensions of the density maps
pixelsperdegree = 642.7000;
paraFovealExtent = 2;

% Create a map that will be used to filter "extreme" values
dParams = [1477 -0.3396 7846 -1.3049 629];
density = @(x) dParams(1).*exp(dParams(2).*x)+dParams(3).*exp(dParams(4).*x)+dParams(5);
[R,C] = ndgrid(1:newDim, 1:newDim);
R = R - newDim/2;
C = C - newDim/2;
r = sqrt(R.^2+C.^2) ./ pixelsperdegree;
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
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
warning('off','signal:findpeaks:largeMinPeakHeight');

% We are going to loop through the processing twice, once each for the
% confocal and split datasets
for dd = 1:2
    
    switch dd
        case 1
            resultFiles=dir('**/confocal/Results_foveated/*_confocal_Fouriest_Result.mat');
            dataFileName = 'confocalDensityProfileData.mat';
            tagName = '_confocal';
            maxThresh = density(r).*3;
        case 2
            resultFiles=dir('**/split detection/Results_foveated/*_split_Fouriest_Result.mat');
            dataFileName = 'splitDensityProfileData.mat';
            tagName = '_split';
            maxThresh = density(r).*2;
    end
    
    % Loop over the result files
    for rr = 1:length(resultFiles)
        
        % Report our progress
        fprintf([resultFiles(rr).name '\n']);
        
        % Load the next file
        fileName = fullfile(resultFiles(rr).folder,resultFiles(rr).name);
        load(fileName);       
        
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
        
        % Extract the subject name
        tmp = strsplit(fileName,filesep);
        subName = tmp{end-3};
        
        % Save the fovea_coords from the confocal, and apply these to the
        % split
        switch dd
            case 1
                % Special case 11061_OD
                if strcmp(subName,'11061_OD')
                    fovea_coords = [4.6105e3 4.75559e3];
                end
                % Special case 11099_OD
                if strcmp(subName,'11099_OD')
                    fovea_coords = [8.5285e3, 7.8643e3];
                end
                foveaCoordStore.(['s_' subName])=fovea_coords;
            case 2
                fovea_coords=foveaCoordStore.(['s_' subName]);
        end        

        % Detect the special case of subject 11051, who wore a -8.5 D
        % spectacle lens during collection of their adaptive optics images.
        % Need to adjust for spectaacle magnification
        if strcmp(subName,'11051_OS') || strcmp(subName,'11051_OD')
            % Need to resize the density_map, and the foveamask, and adjust
            % the fovea_coords
            sg = createSceneGeometry(...
                'sphericalAmetropia',-8.5000,...
                'axialLength',25.9250,...
                'spectacleLens',-8.5);
            magFactor = 1/sg.refraction.cameraToRetina.magnification.spectacle;
            density_map = imresize(density_map,magFactor);
            foveamask = round(imresize(foveamask,magFactor));
            fovea_coords = fovea_coords.*magFactor;
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
        figHandle = figure('Name',subName,'Position',  [100, 100, 800, 200]);
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

        % Show the filtered map
        subplot(1,3,2)
        imagesc(imDensity)
        hold on
        plot([1 newDim],round([newDim newDim]/2),'-r');
        plot(round([newDim newDim]/2),[1 newDim],'-r');
        axis square
        axis off
        title('filter periphery')
        
        % Set up the fovea mask
        imFovea = zeros(newDim,newDim);
        imFovea(1:imsize(1),1:imsize(2))=single(~foveamask);
        imFovea=imtranslate(imFovea,offset);
        
        % Dilate the foveamask
        imFovea = imdilate(imFovea,strel('square',foveaDilate));
        
        % Down-sample the maps
        imDensity = imresize(imDensity,downSample);
        imFovea = imresize(imFovea,downSample);
        
        % Flip left eye density maps so they are pseudo-right eye
        if strcmp(laterality,'OS')
            imDensity = fliplr(imDensity);
            imFovea = fliplr(imFovea);
        end
        
        % Create polar maps
        polarDensity = convertImageMapToPolarMap(imDensity);
        polarFovea = convertImageMapToPolarMap(imFovea);
        polarDim = newDim*downSample*2-1;

        % Show the polar density map pre filtering
        subplot(1,3,3)
        imagesc(polarDensity)
        axis square
        axis off
        title('filter center; polar')

        % Calculate the support in degrees for the polar image
        supportDeg = (1:polarDim)./(pixelsperdegree*downSample*4);        
        
        % Remove decreasing components close to the fovea in the confocal
        % images
        supportDegIdx = find(supportDeg>paraFovealExtent,1);
        threshIdx = ones(1,size(polarDensity,1));
        
        if dd==1
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
                
            end
        end
        
        % Add the filtering ridge
        hold on
        plot(threshIdx,1:size(polarDensity,1),'-r');
        drawnow
        
        % Apply the filtering
        for nn=1:size(polarDensity,1)
            polarDensity(nn,1:threshIdx(nn))=nan;
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
        data.meta.foveaDilate = foveaDilate;
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

% Turn on a warning
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

