dirRootPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/2DThicknessMapsAllLayers_10_19_2018';
subjectList = dir(dirRootPath);
subjectList = subjectList(4:end);

nSubs = length(subjectList);

layerIdx = [2,3];
thickness=[];
allSubs = [];
imageSize = [768 768];
foveaThickThresh = 75;
everyThicknessMap = nan(nSubs,2,imageSize(1),imageSize(2));

for ss=1:nSubs
    
    % Assemble the path to the data directory for this subject
    dirPath = fullfile(subjectList(ss).folder,subjectList(ss).name);
    
    % Obtain the list of result files for this subject
    fileList = dir(fullfile(dirPath,'*.mat'));
    
    % If there are no files for this subject, move on to the next subject
    if isempty(fileList)
        continue
    end
    
    % If there are more than two result files for this subject, error.
    if length(fileList)>2
        error('Why are more than two eyes for this subject?')
    end
    
    % Report this subject's name
    fprintf(['Processing subject: ' subjectList(ss).name '\n']);
    
    % Open a figure to display the maps
    %    figure('Name',subjectList(ss).name,'NumberTitle','off');
    
    % Loop over the right and left eye
    for ii = 1:length(fileList)
        
        % There are a few eyes with bad registrations. We skip them
        % here until such time as Min has them working again.
        if any(contains(fileList(ii).name,{'11028_OS.mat','11055_OS.mat','11088_OD.mat','11091_OS.mat'}))
            continue
        end
        
        % Assemble the name of this file and load it
        fileName = fullfile(fileList(ii).folder,fileList(ii).name);
        load(fileName);
        
        % Get the thickness by summing the layers defined in
        % layerIdx
        thisThick = sum(subject.BothMeanLayerThicknessesOnSLOInterp(:,:,layerIdx),3);
        
        % If these are data from the left eye, mirror reverse
        if contains(fileList(ii).name,'_OS.mat')
            thisThick = fliplr(thisThick);
        end
        
        % Set points with zero thickness to nan
        thisThick(thisThick==0)=nan;
        
        % Resize the map to ensure that all maps have the same
        % dimensions (some were acquired in 1536x1536 resolution
        thisThick=imresize(thisThick,imageSize,'bilinear');
        
        % Now shift the image so that the fovea is at a standard
        % location. First, find the fovea at the center of the
        % image. This is defined by taking the weighted mean of the
        % thinnest portion of the central 20% of the image
        w=thisThick;
        w(1:round(imageSize(1)*0.4),:)=nan;
        w(round(imageSize(1)*0.6):end,:)=nan;
        w(:,1:round(imageSize(1)*0.4))=nan;
        w(:,round(imageSize(1)*0.6):end)=nan;
        w(w>foveaThickThresh)=nan;
        w=max(max(w))-w;
        w=w(1:prod(imageSize));
        w=w./nansum(w);
        [X,Y] = ind2sub(imageSize,1:prod(imageSize));
        foveaCoord(1) = nansum(X.*w(1:prod(imageSize)));
        foveaCoord(2) = nansum(Y.*w(1:prod(imageSize)));
        
        % Save a mask of the nan values
        nanMask = zeros(size(thisThick));
        nanMask(isnan(thisThick))=1;
        thisThick(isnan(thisThick))=0;
        
        % Shift the thickness and nanmask images to align the foveas
        thisThick = imtranslate(thisThick,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
        nanMask = imtranslate(nanMask,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
        
        % Restore the nans
        thisThick(nanMask>0.25)=nan;
        
        % Store the this thickness map in the full array
        if contains(fileList(ii).name,'_OS.mat')
            everyThicknessMap(ss,1,:,:)=thisThick;
        else
            everyThicknessMap(ss,2,:,:)=thisThick;
        end
        
        % Show this map
        %             subplot(1,2,ii)
        %             imagesc(thisThick);
        %             axis square
        %             drawnow
        
    end % Looping over the two eyes
end % Looping over subjects

% Create the grand average map. First take the nanmean within
% subject across eyes. This takes the union of all available measurements.
% Then, take the mean across subjects. This finds the intersection across
% subjects of the locations with data.
grandAverage = squeeze(mean(nanmean(everyThicknessMap,2),1));

% Trim the optic disc area away as this region contains segmentation
% artifacts
grandAverage(:,640:end)=nan;

% Display the grand average
figure
mesh(grandAverage)
xlim([0 imageSize(1)]);
ylim([0 imageSize(1)]);
axis square

% Find the set of map locations for which every subject has a measurement
observedIdx = ~isnan(grandAverage);

% Obtain the mean thickness for the left and right eye for each subject,
% plot and report the across subject correlation of the two eyes
for ss=1:nSubs
    for ii=1:2
        tmp = squeeze(everyThicknessMap(ss,ii,:,:));
        meanThickByEye(ss,ii) = nanmean(tmp(observedIdx));
    end
end
figure
plot(meanThickByEye(:,1),meanThickByEye(:,2),'*r');
refline(1,0);
xlim([45 70])
ylim([45 70])
axis square
xlabel('RGC+IPL mean thickness OS');
ylabel('RGC+IPL mean thickness OD');
rho = corrcoef(meanThickByEye(:,1),meanThickByEye(:,2),'rows','pairwise');
fprintf('Left and right RGC+IPL thickness correlation across subjects R (%d df) = %2.0f \n',sum(~isnan(mean(meanThickByEye,2))),rho(1,2));

% Conduct a PCA analysis across subjects, averaged over eyes. First,
% Convert the non-nan portion of the image to a vector
X = [];
avgMapBySubject = squeeze(nanmean(everyThicknessMap,2));
for ss=1:nSubs
    tmp = squeeze(avgMapBySubject(ss,:,:));
    X(ss,:) = tmp(observedIdx);
end

% Calc the PCA
[coeff,score,latent,tsquared,explained,mu] = pca(X,'Centered',false);
Xfit = score(:,1:3)*coeff(:,1:3)';

% Create a 3D plot of the scores for the first 3 components
figure
plot3(score(:,1),score(:,2),score(:,3),'*r')

% Now show images of the first, second, and third PCA components
coeffMap = nan(imageSize(1),imageSize(2));
figure
for ii=1:4
    subplot(2,2,ii);
    coeffMap(observedIdx)=coeff(:,ii);
    imagesc(coeffMap);
    xlim([0 imageSize(1)]);
    ylim([0 imageSize(1)]);
    axis square
end


% mesh(thickness)
%
% dims = size(thickness);
% [x,y] = ind2sub(dims,1:prod(dims));
% z = thickness(1:prod(dims));
% gf = gridfit(x,y,z,0:10:dims(1),0:10:dims(2),'smoothness',1,'interp','bilinear');
% surf(0:10:dims(1),0:10:dims(2),gf);
%
%
