dirRootPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/2DThicknessMapsAllLayers_10_19_2018';
subjectList = dir(dirRootPath);

layerIdx = [2,3];
thickness=[];
allSubs = [];
imageSize = [768 768];
foveaThickThresh = 75;

for ss=3:length(subjectList)
    
    dirPath = fullfile(subjectList(ss).folder,subjectList(ss).name);
    
    fileList = dir(fullfile(dirPath,'*.mat'));
    if ~isempty(fileList)
        subjectList(ss).name
        thisThick = [];
        allThick = [];
        for ii = 1:length(fileList)
            fileName = fullfile(fileList(ii).folder,fileList(ii).name);
            load(fileName)
            thisThick = sum(subject.BothMeanLayerThicknessesOnSLOInterp(:,:,layerIdx),3);
            if contains(fileList(ii).name,'_OS.mat')
                thisThick = fliplr(thisThick);
            end
            allThick(ii,:,:)=thisThick;
        end
        
        allThick(allThick==0)=nan;
        thickness = squeeze(mean(allThick,1));
        thickness=imresize(thickness,imageSize,'bilinear');
        
        
        % Find the fovea at the center of the image
        w=thickness;
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
        nanMask = zeros(size(thickness));
        nanMask(isnan(thickness))=1;
        thickness(isnan(thickness))=0;
        
        % Shift the thickness and nanmask images to align the foveas
        thickness = imtranslate(thickness,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
        nanMask = imtranslate(nanMask,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
        
        % Restore the nans
        thickness(nanMask>0.25)=nan;
        
        % Store the thickness map for this subject
        allSubs(ss,:,:)=thickness;
    end
end


% Remove empty entries from allSubs
allSubs(allSubs==0)=nan;
allNanEntries = squeeze(all(all(isnan(allSubs),2),3));
allSubs = allSubs(~allNanEntries,:,:);

% Trim the optic disc area away
allSubs(:,:,640:end)=nan;

% Make the average and SD map
final = squeeze(nanmean(allSubs,1));
finalSD = squeeze(nanstd(allSubs,1));

figure
mesh(final)
figure
mesh(finalSD)

% Conduct a probabilistic PCA analysis
% Downsample and vectorize the RGC thickness maps
for ii=1:size(allSubs,1)
    X(ii,:) = reshape( imresize(squeeze(allSubs(ii,:,:)),round(imageSize./20),'bilinear'), [1 prod(round(imageSize./20))] );
end

% Remove any columns that are all nans
nanCols=all(isnan(X),1);
Xsub = X(:,~nanCols);
options=statset();
options.Display='iter';
[coeffSub,score,pcvar,mu] = ppca(Xsub,3,'Options',options);
coeff(~nanCols,~nanCols)=coeffSub;


% mesh(thickness)
%
% dims = size(thickness);
% [x,y] = ind2sub(dims,1:prod(dims));
% z = thickness(1:prod(dims));
% gf = gridfit(x,y,z,0:10:dims(1),0:10:dims(2),'smoothness',1,'interp','bilinear');
% surf(0:10:dims(1),0:10:dims(2),gf);


