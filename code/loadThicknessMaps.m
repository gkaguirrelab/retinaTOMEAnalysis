dirRootPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/2DThicknessMapsAllLayers_10_19_2018';
rgcMapSavePath = '/Users/aguirre/Documents/MATLAB/projects/rgcPopulationModel/data/rgcIplThicknessMap.mat';

subjectList = dir(dirRootPath);
subjectList = subjectList(4:end);

nSubs = length(subjectList);

layerSetLabels = {'RGC+IPL','RNFL','OPL','total'};
layerSets = [{2:3},{1},{6},{1:11}];
nSets = length(layerSets);
rgcIPLidx = 1;
layerIdx = [2,3];
thickness=[];
allSubs = [];
imageSize = [768 768];
foveaThickThresh = 75;
everyThicknessMap = nan(nSets,nSubs,2,imageSize(1),imageSize(2));


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
    figure('Name',subjectList(ss).name,'NumberTitle','off');
    
    if length(fileList)==1
        warning(['subject ' subjectList(ss).name 'has just one eye']);
    end
    
    % Loop over the right and left eye
    for ii = 1:length(fileList)
                
        % Assemble the name of this file and load it
        fileName = fullfile(fileList(ii).folder,fileList(ii).name);
        load(fileName);
        
        % Empty the variable that will hold the fovea coordinates
        foveaCoord=[];
        
        % Loop over the layer sets
        for jj=1:nSets
            
            % Get the thickness by summing the layers defined in
            % layerIdx
            thisThick = sum(subject.BothMeanLayerThicknessesOnSLOInterp(:,:,layerSets{jj}),3);
            
            % If these are data from the left eye, mirror reverse
            if contains(fileList(ii).name,'_OS.mat')
                thisThick = fliplr(thisThick);
            end
            
            % Set points with zero thickness to nan
            thisThick(thisThick==0)=nan;
            
            % Resize the map to ensure that all maps have the same
            % dimensions (some were acquired in 1536x1536 resolution
            thisThick=imresize(thisThick,imageSize,'bilinear');
            
            % If this is the first layer set, then we are working with the
            % RGC-IPL layer. Use this to find the fovea at the center of the
            % image. This is defined by taking the weighted mean of the
            % thinnest portion of the central 20% of the image
            if jj==1
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
            end
            
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
                everyThicknessMap(jj,ss,1,:,:)=thisThick;
            else
                everyThicknessMap(jj,ss,2,:,:)=thisThick;
            end
            
            % Show the RGC+IPL map
            if jj==1
                tmp = fliplr(thisThick);
                if contains(fileList(ii).name,'_OS.mat')
                    subplot(1,2,1)
                    tmp = fliplr(tmp);
                else
                    subplot(1,2,2)
                end
                imagesc(tmp);
                hold on
                plot(imageSize(1)/2,imageSize(2)/2,'+k')
                axis square
                drawnow
            end
            
        end % Loop over layer sets
        
    end % Looping over the two eyes
end % Looping over subjects

% Create the grand averages map for each set. First take the nanmean within
% subject across eyes. This takes the union of all available measurements.
% Then, take the mean across subjects. This finds the intersection across
% subjects of the locations with data.
grandAverage = squeeze(mean(nanmean(everyThicknessMap,3),2));

% Trim away areas with artifacts and strangeness at the edge of the image.
% The largest of these is in the area of the optic disc
grandAverage(:,:,640:end)=nan;
grandAverage(:,:,640:end)=nan;
grandAverage(:,1:20,:)=nan;
grandAverage(:,1:20,:)=nan;
grandAverage(:,712:end,:)=nan;

% Display the grand averages
figure
for jj=1:nSets
    subplot(2,2,jj)
    mesh(squeeze(grandAverage(jj,:,:)));
    xlim([0 imageSize(1)]);
    ylim([0 imageSize(1)]);
    axis square
end

% Save the RGC+IPL thickness map
rgcIplThicknessMap = squeeze(grandAverage(1,:,:));
save(rgcMapSavePath,'rgcIplThicknessMap');

% Find the set of map locations for which every subject has a measurement
observedIdx = ~isnan(squeeze(sum(grandAverage,1)));

% Obtain the mean thickness for the left and right eye for each subject,
for ss=1:nSubs
    for ii=1:2
        for jj=1:nSets
            tmp = squeeze(everyThicknessMap(jj,ss,ii,:,:));
            meanThickByEye(jj,ss,ii) = nanmean(tmp(observedIdx));
        end
    end
end

% Plot and report the across subject correlation of the two eyes
figure
for jj=1:nSets
    subplot(2,2,jj);
    plot(meanThickByEye(jj,:,1),meanThickByEye(jj,:,2),'*r');
    refline(1,0);
    axis equal
    axis square
    xlabel([layerSetLabels{jj} ' mean thickness OS']);
    ylabel([layerSetLabels{jj} ' mean thickness OD']);
    rho = corrcoef(meanThickByEye(jj,:,1),meanThickByEye(jj,:,2),'rows','pairwise');
    fprintf(['Left and right ' layerSetLabels{jj} ' thickness correlation across subjects R (%d df) = %2.2f \n'],sum(~isnan(mean(squeeze(meanThickByEye(jj,:,:)),2))),rho(1,2));
end

% Report the correlation coefficients of the measures (averaged across
% eyes) across subjects
layerSetLabels
corrcoef(nanmean(meanThickByEye,3)','rows','pairwise')

% Conduct a PCA analysis across subjects, averaged over eyes. First,
% Convert the non-nan portion of the image to a vector
X = [];
avgMapBySubject = squeeze(nanmean(squeeze(everyThicknessMap(rgcIPLidx,:,:,:,:)),2));
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
    mesh(coeffMap);
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
