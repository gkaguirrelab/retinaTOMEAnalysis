clear
%close all

downSample = 0.1;
foveaDilate = 100;
angleAccumulate = 45;
armFilterPointDegrees = 4;
armFilterDensityThresh = 2000;
recalculateFlag = false;


% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')

% Find the list of result files
%resultFiles=dir('**/confocal/Results/*_confocal_Fouriest_Result.mat');
%dataFileName = 'densityProfileData.mat';

resultFiles=dir('**/split detection/Results/*_split_Fouriest_Result.mat');
dataFileName = 'splitDensityProfileData.mat';


% Either load or create the data variable
if ~recalculateFlag
    
    load(dataFileName)
    
else
    
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
            filterIdx = find(supportDeg > armFilterPointDegrees,1);
            filterRegion = dataRectangle(:,filterIdx:end);
            filterRegion(filterRegion(:)>armFilterDensityThresh)=nan;
            dataRectangle(:,filterIdx:end) = filterRegion;
            meridianDensity(mm,:)=nanmean(dataRectangle);
            
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


%% Aggregate the profiles
% Find the longest support deg
[supportLength,idx] = max(cellfun(@(x) length(x.profile.supportDeg),data));
supportDeg = 0:data{1}.meta.supportDegDelta:data{1}.meta.supportDegDelta*(supportLength-1);
polarToMeridian = {'Inferior','Nasal','Superior','Temporal'};

badIdx = zeros(size(data));
data = data(logical(~badIdx));

% Loop over the arms
for mm = 1:length(polarToMeridian)
    dataMatrix = nan(length(data),supportLength);
    for rr = 1:length(data)
        tmp = data{rr}.profile.(polarToMeridian{mm});
        dataMatrix(rr,1:length(tmp))=tmp;
    end
    dataAggregate.(polarToMeridian{mm}).median = nanmedian(dataMatrix);
    dataAggregate.(polarToMeridian{mm}).stdev = nanstd(dataMatrix);
end
dataAggregate.supportDeg = supportDeg;


%% Display profiles
figure
startIdx = find(supportDeg > 0.5,1);
for mm = 1:length(polarToMeridian)
    plot(dataAggregate.supportDeg(startIdx:end),dataAggregate.(polarToMeridian{mm}).median(startIdx:end));
    hold on
end
title(sprintf('Across-subject median density in %d degree wedge',angleAccumulate))
ylim([0 7000]);
xlim([0 15]);
ylabel('mean density [cones/deg^2]');
xlabel('distance from fovea [deg]')
legend(polarToMeridian,'FontSize',16)




%
%
% idxToDisplay = 3;
% imDensity = data{idxToDisplay}.imDensity;
% polarDensity = data{idxToDisplay}.polarDensity;
%
% %% Display maps
% figure
% tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
%
% nexttile
% imagesc(imDensity)
% axis equal
% axis off
% title('density map')
%
% nexttile
% imagesc(polarDensity)
% axis equal
% axis off
% title('density polar angle')
% hold on
% for mm = 1:4
%     if mm==4
%         plot([1 1],[0 meridianWidth],'-r','LineWidth',2);
%         plot([1 1],[polarDim-meridianWidth polarDim],'-r','LineWidth',2);
%         plot([1 polarDim],[meridianWidth meridianWidth],'-r');
%         plot([1 polarDim],[polarDim-meridianWidth polarDim-meridianWidth],'-r');
%     else
%         plot([1 1],[meridianIdx*mm-meridianWidth meridianIdx*mm+meridianWidth],'-r','LineWidth',2);
%         plot([1 polarDim],[meridianIdx*mm-meridianWidth meridianIdx*mm-meridianWidth],'-r');
%         plot([1 polarDim],[meridianIdx*mm+meridianWidth meridianIdx*mm+meridianWidth],'-r');
%     end
% end

%

