%% set input parameters
radMM=3;
smpPerMM=16;
sectorAngle=3;

%% Load dispMap
dispMap = makeMap(radMM,smpPerMM,sectorAngle);

%% loop over rotations

rotDegs =0:sectorAngle:360-sectorAngle;
sampleBase = 0:1/smpPerMM:radMM;
for i = 1:length(rotDegs)
    %% rotate and crop
    disp_crop   = imrotate(dispMap,-1.*rotDegs(i),'crop','nearest');
    %% Get meridian from RF density rotated input image
    disp_midPoint = round(size(disp_crop,1)/2);
    dispVec(i,:) = disp_crop(disp_midPoint,disp_midPoint:end);
    newLocVec(i,:) = sampleBase - dispVec(i,:);
    diffVec(i,:) = diff(newLocVec(i,:));
end






% 
% sampleBase = 0:1/smpPerMM:radMM;
% sampleBaseMat = repmat(sampleBase,[size(dispVec,1),1]);
% 
% newLocation = sampleBaseMat - dispVec;
% 
% 
% figure;
% subplot(1,3,1) 
% plot(0:1/16:3,dispVec(1,:))
% title('Displacement vectors along the nasal')
% subplot(1,3,2) 
% plot(0:1/16:3,newLocation(1,:))
% title('New location of along nasal at ecc')
% subplot(1,3,3) 
% diffVec= diff(newLocation,2);
% plot(0:1/16:3,diffVec(1,:))
% ylim([-.1 .1])
% 
% title('Diff of new location')
% 
% 
% 
% figure;
% subplot(1,3,1) 
% imagesc(dispVec);
% title('Displacement vectors along the radial')
% subplot(1,3,2) 
% imagesc(newLocation)
% 
% title('New location of index at an ecc')
% 
% subplot(1,3,3) 
% imagesc(diff(newLocation))
% 
% title('Diff of new location')