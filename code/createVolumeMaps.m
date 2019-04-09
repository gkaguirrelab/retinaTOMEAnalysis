inDirThickness = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\averageThicknessMapsBySubject';
inDirmmPerDeg = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\mmPerDegMaps';
outDir = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\volumeMapsBySubject';

%find all subjects
subIDs = dir(fullfile(inDirThickness,'1*'));

for i = 1:length(subIDs)

%thicknessMap = squeeze(mean(everyThicknessMap(1,5,:,:,:),3));

%load mmPerDeg Maps
LoadmmPerDeg=load(fullfile(inDirmmPerDeg,[subIDs(i).name '_mmPerDegMap.mat']));
mmPerDeg=LoadmmPerDeg.mmPerDeg;
%load average thickness maps
LoadthicknessMap=load(fullfile(inDirThickness,subIDs(i).name,[subIDs(i).name '_averageMaps.mat']));
%LoadthicknessMap.averageMaps;


mkdir(fullfile(outDir, subIDs(i).name));

for L = 1:4
    switch L
        case 1
            thicknessMap = LoadthicknessMap.averageMaps.RGCIPL;
            savename = fullfile(outDir, subIDs(i).name, [subIDs(i).name '_RGCIPL_volumeMaps.mat']);
        case 2
            thicknessMap = LoadthicknessMap.averageMaps.RNFL;
            savename = fullfile(outDir, subIDs(i).name, [subIDs(i).name '_RNFL_volumeMaps.mat']);
        case 3
            thicknessMap = LoadthicknessMap.averageMaps.OPL;
            savename = fullfile(outDir, subIDs(i).name, [subIDs(i).name '_OPL_volumeMaps.mat']);
        case 4
            thicknessMap = LoadthicknessMap.averageMaps.TotalRetina;
            savename = fullfile(outDir, subIDs(i).name, [subIDs(i).name '_TotalRetina_volumeMaps.mat']);
    end
    
XN = size(thicknessMap,1);
YN = size(thicknessMap,2);

%fill in nans
mmPerDegInterp = fillmissing(mmPerDeg,'linear');

%roorient image to OD clinical view
mmPerDegInterp_rot = fliplr(mmPerDegInterp');

%resize to same size as slo
mmPerDegMapInterp_rot_resize = imresize(mmPerDegInterp_rot,[XN YN]);

%find size of pixel in deg^2
%30deg by 30deg fov divided by XN pixxels x YN pixels
degsquredperpixel = (30/XN)*(30/YN);

%convert thickness map to volume map with thickness(mm)*degree^2
volumeMap_mmDegSquared = thicknessMap*degsquredperpixel;

%convert volumes from thickness(mm)*degree^2 to mm^3
volumeMap_mmCubed = volumeMap_mmDegSquared.*mmPerDegMapInterp_rot_resize.*mmPerDegMapInterp_rot_resize;

save(savename,'mmPerDegMapInterp_rot_resize','volumeMap_mmDegSquared','volumeMap_mmCubed');

end

% figure(1)
% imshow(mmPerDegMapInterp_rot_resize)
% caxis([min(mmPerDegMapInterp_rot_resize(:)) max(mmPerDegMapInterp_rot_resize(:))])
% colorbar
% figure(2)
% imshow(volumeMap_mmDegSquared)
% caxis([min(volumeMap_mmDegSquared(:)) max(volumeMap_mmDegSquared(:))])
% colorbar
% figure(3)
% imshow(volumeMap_mmCubed) 
% caxis([min(volumeMap_mmCubed(:)) max(volumeMap_mmCubed(:))])
% colorbar

end
