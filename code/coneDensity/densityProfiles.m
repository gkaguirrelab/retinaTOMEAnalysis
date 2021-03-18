cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/11015_OD/confocal/Results')

load('NC_11015_20160401_OD_confocal_Fouriest_Result.mat');


downSample = 0.1;
foveaDilate = 100;
angleAccumulate = 45;

pixelsperdegree;
fovea_coords;
density_map;

% Pad the density map so that it is square around the fovea_coords
imsize = size(density_map);
newDim =  max(max([imsize-round(fovea_coords);round(fovea_coords)]))*2;
offset = round(newDim/2-fovea_coords);

% Set up the density map
imDensity = zeros(newDim,newDim);
imDensity(1:imsize(1),1:imsize(2))=density_map;
imDensity(isnan(imDensity))=0;
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

% Mask the polarDensity by the polarFovea
polarDensity(polarFovea > 0.1) = nan;

% Calculate the support in degrees for the polar image
supportDeg = (1:newDim*downSample*2)./(pixelsperdegree*downSample*2);
supportDeg = supportDeg(2:end);

% Loop through the merdians
merdianDensity = nan(4,newDim*downSample*2-1);
meridianIdx = floor(newDim*downSample*2/4);
meridianWidth = round(((newDim*downSample*2-1)/360)*angleAccumulate/2);
for mm = 1:4
    if mm==4
        merdianDensity(mm,:)=nanmean( [ polarDensity(1:meridianWidth,:); ...
            polarDensity(end-meridianWidth:end,:)] );
    else
        merdianDensity(mm,:)=nanmean(polarDensity(meridianIdx*mm-meridianWidth:meridianIdx*mm+meridianWidth,:));
    end
end

%% Display maps
figure
tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');

nexttile
imagesc(imDensity)
axis equal
axis off
title('density map')

nexttile
imagesc(imFovea)
axis equal
axis off
title('fovea map')

nexttile
imagesc(polarDensity)
axis equal
axis off
title('density polar angle')
hold on
for mm = 1:4
    if mm==4
        plot([1 1],[0 meridianWidth],'-r','LineWidth',2);
        plot([1 1],[newDim*downSample*2-1-meridianWidth newDim*downSample*2-1],'-r','LineWidth',2);        
    else
        plot([1 1],[meridianIdx*mm-meridianWidth meridianIdx*mm+meridianWidth],'-r','LineWidth',2);
    end
end

nexttile
imagesc(polarFovea)
axis equal
axis off
title('fovea polar angle')


%% Display profiles
figure
plot(supportDeg,merdianDensity);
title(sprintf('Density in %d degree wedge',angleAccumulate))
ylabel('mean density [cones/deg^2]');
xlabel('distance from fovea [deg]')
legend({'0°','90°','180°','270°'},'FontSize',24)
