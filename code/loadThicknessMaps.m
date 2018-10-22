dirPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/2DThicknessMapsAllLayers_10_19_2018/11015';
fileList = dir(fullfile(dirPath,'*.mat'));

layerIdx = [2,3];

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
thickness = squeeze(sum(allThick,1));
mesh(thickness)

dims = size(thickness);
[x,y] = ind2sub(dims,1:prod(dims));
z = thickness(1:prod(dims));
gf = gridfit(x,y,z,0:10:dims(1),0:10:dims(2),'smoothness',1,'interp','bilinear');