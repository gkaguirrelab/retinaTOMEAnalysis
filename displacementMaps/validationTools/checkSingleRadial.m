function checkSingleRadial(map,radialRotDeg,sizeMM,smpPerMM)

%Create sample base in MM based on map size and sample rate
sampleBase=0:1/smpPerMM:sizeMM;

rotMap   = imrotate(map,-1.*radialRotDeg,'crop','nearest');

center = round(size(rotMap,1));

radial = rotMap(center,center:end);

newPos = sampleBase - radial;

diffRadial = diff(

figure;
subplot(
