clear

% Hard coded values
downSample = 0.1;
foveaDilate = 100;
angleAccumulate = 45;
armFilterPointsDegrees = [0.5, 4];
armFilterDensityThresh = [10000, 2000];

% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')

% We are going to loop through the processing twice, once each for the
% confocal and split datasets
for dd = 1:2
    
    switch dd
        case 1
            resultFiles=dir('**/confocal/Results/*_confocal_Fouriest_Result.mat');
            dataFileName = 'confocalDensityProfileData.mat';
        case 2
            resultFiles=dir('**/split detection/Results/*_split_Fouriest_Result.mat');
            dataFileName = 'splitDensityProfileData.mat';
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
        
        % Pad the density map so that it is square around the fovea_coords
        imsize = size(density_map);
        newDim =  max(max([imsize-round(fovea_coords);round(fovea_coords)]))*2;
        newDim = round(newDim+9,-1); % Ensure that it is a multiple of 10
        offset = round(newDim/2-fovea_coords);
        
        % Set up the density map
        imDensity = zeros(newDim,newDim);
        imDensity(1:imsize(1),1:imsize(2))=density_map;
        imDensity(isnan(imDensity))=0;
        imDensity(isinf(imDensity))=0;
        imDensity=imtranslate(imDensity,offset);
        imDensity(imDensity==0)=nan;
        
        % Set up the fovea mask
        imFovea = zeros(newDim,newDim);
        imFovea(1:imsize(1),1:imsize(2))=single(~foveamask);
        imFovea=imtranslate(imFovea,offset);
        
        % Dilate the foveamask
        imFovea = imdilate(imFovea,strel('square',foveaDilate));
        
        % Down-sample the maps
        imDensity = imresize(imDensity,downSample);
        imFovea = imresize(imFovea,downSample);
        
        % Create polar maps
        polarDensity = convertImageMapToPolarMap(imDensity);
        polarFovea = convertImageMapToPolarMap(imFovea);
        polarDim = newDim*downSample*2-1;
        
        % Mask the polarDensity by the polarFovea
        polarDensity(polarFovea > 0.1) = nan;
        
        % Calculate the support in degrees for the polar image
        supportDeg = (1:polarDim)./(pixelsperdegree*downSample*4);
        
        % Loop through the merdians and grab the data rectangle
        meridianDensity = nan(4,polarDim);
        meridianIdx = floor(polarDim/4);
        meridianWidth = round(((polarDim)/360)*angleAccumulate/2);
        for mm = 1:4
            if mm==4
                dataRectangle=[ polarDensity(1:meridianWidth,:); ...
                    polarDensity(end-meridianWidth:end,:)];
            else
                dataRectangle = polarDensity(meridianIdx*mm-meridianWidth:meridianIdx*mm+meridianWidth,:);
            end
            
            % Filter out unreasonable values from beyond the armFilterPointDegrees
            for ff = 1:length(armFilterPointsDegrees)
                filterIdx = find(supportDeg > armFilterPointsDegrees(ff),1);
                filterRegion = dataRectangle(:,filterIdx:end);
                filterRegion(filterRegion(:)>armFilterDensityThresh(ff))=nan;
                dataRectangle(:,filterIdx:end) = filterRegion;
                meridianDensity(mm,:)=nanmean(dataRectangle);
            end
        end
        
        % Store the data
        data{rr}.meta.filename = resultFiles(rr).name;
        data{rr}.meta.folder = resultFiles(rr).folder;
        data{rr}.meta.laterality = laterality;
        data{rr}.meta.downSample = downSample;
        data{rr}.meta.foveaDilate = foveaDilate;
        data{rr}.meta.angleAccumulate = angleAccumulate;
        data{rr}.meta.supportDegDelta = supportDeg(1);
        data{rr}.imDensity = imDensity;
        data{rr}.polarDensity = polarDensity;
        data{rr}.profile.supportDeg = supportDeg;
        for mm = 1:4
            data{rr}.profile.(polarToMeridian{mm}) = meridianDensity(mm,:);
        end
        
    end
    
    % Save the assembled data file
    save(dataFileName,'data','-v7.3');
    
end
