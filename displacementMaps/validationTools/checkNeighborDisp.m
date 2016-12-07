%% set input parameters
radMM=3;
smpPerMM=16;
sectorAngle=3;

%% Load dispMap
dispMap = makeMap(radMM,smpPerMM,sectorAngle);

%% loop over rotations

rotDegs =0:sectorAngle:360-sectorAngle;

for i = 1:length(rotDegs)
    %% rotate and crop
    disp_crop   = imrotate(dispMap,-1.*rotDegs(i),'crop','bilinear');
    %% Get meridian from RF density rotated input image
    disp_midPoint = round(size(disp_crop,1)/2);
    dispVec(i,:) = disp_crop(disp_midPoint,disp_midPoint:end);
end

sampleBase = 0:1/smpPerMM:radMM;
sampleBaseMat = repmat(sampleBase,[size(dispVec,1),1]);

newLocation = sampleBaseMat - dispVec;
