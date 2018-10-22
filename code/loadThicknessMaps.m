dirRootPath = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis/2DThicknessMapsAllLayers_10_19_2018';
subjectList = dir(dirRootPath);

layerIdx = [2,3];
thickness=[];
allSubs = [];
imageSize = [768 768];

for ss=3:length(subjectList)
    
    dirPath = fullfile(subjectList(ss).folder,subjectList(ss).name);
    
    fileList = dir(fullfile(dirPath,'*.mat'));
    if ~isempty(fileList)
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
        
        allSubs(ss,:,:)=thickness;
        figure
        imagesc(thickness);
        subjectList(ss).name
        drawnow
    end
end



% mesh(thickness)
%
% dims = size(thickness);
% [x,y] = ind2sub(dims,1:prod(dims));
% z = thickness(1:prod(dims));
% gf = gridfit(x,y,z,0:10:dims(1),0:10:dims(2),'smoothness',1,'interp','bilinear');
% surf(0:10:dims(1),0:10:dims(2),gf);