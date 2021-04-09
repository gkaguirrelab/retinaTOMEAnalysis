clear

% Hard coded values
downSample = 0.05;
foveaDilate = 100;
newDim = 18000; % Dimensions of the density maps
pixelsperdegree = 642.7000;

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

cd(fullfile(dropboxBaseDir,'Connectome_AOmontages_images'))

% Turn off a warning
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');

% We are going to loop through the processing twice, once each for the
% confocal and split datasets
for dd = 1:2
    
    switch dd
        case 1
            resultFiles=dir('**/confocal/Results_foveated/*_confocal_Fouriest_Result.mat');
            dataFileName = 'confocalDensityProfileData.mat';
            maxThresh = density(r).*3;
        case 2
            resultFiles=dir('**/split detection/Results_foveated/*_split_Fouriest_Result.mat');
            dataFileName = 'splitDensityProfileData.mat';
            maxThresh = density(r).*2;
    end
    
    % Create variables to hold the data
    data = cell(1,length(resultFiles));
    
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
        
        % Filter the density map to remove extreme values
        figure('Name',subName,'Position',  [100, 100, 600, 200]);
        subplot(1,3,1)
        imagesc(imDensity)
        axis square
        axis off
        imDensity(imDensity>maxThresh)=nan;
        subplot(1,3,2)
        imagesc(imDensity)
        axis square
        axis off
        
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
        
        % Mask the polarDensity by the polarFovea
        polarDensity(polarFovea > 0.1) = nan;

        % Show the polar density map
        subplot(1,3,3)
        imagesc(polarDensity)
        axis square
        axis off
        drawnow

        % Calculate the support in degrees for the polar image
        supportDeg = (1:polarDim)./(pixelsperdegree*downSample*4);        
        
        % Store the data
        data{rr}.meta.subName = subName;
        data{rr}.meta.filename = resultFiles(rr).name;
        data{rr}.meta.folder = resultFiles(rr).folder;
        data{rr}.meta.laterality = laterality;
        data{rr}.meta.downSample = downSample;
        data{rr}.meta.foveaDilate = foveaDilate;
        data{rr}.meta.supportDegDelta = supportDeg(1);
        data{rr}.imDensity = imDensity;
        data{rr}.polarDensity = polarDensity;
        
    end
    
    % Save the assembled data file
    save(dataFileName,'data','-v7.3');

end

% Turn on a warning
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

