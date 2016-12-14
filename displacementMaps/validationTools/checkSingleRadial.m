function checkSingleRadial(radialRotDeg,sizeMM,smpPerMM,sectorAngle)


map = makeMap(sizeMM,smpPerMM,sectorAngle)
%Create sample base in MM based on map size and sample rate
sampleBase=0:1/smpPerMM:sizeMM;

rotMap   = imrotate(map,-1.*radialRotDeg,'crop','nearest');

center = round(size(rotMap,1)/2);

radial = rotMap(center,center:end);

newPos = sampleBase - radial;

diffRadial = diff(newPos);

figure;
subplot(1,3,1)
plot(sampleBase,radial)
xlabel('Eccentricity (mm)')
ylabel('Displacement (mm)')
title('Displacement Along Radial')
axis square
subplot(1,3,2)
plot(sampleBase,newPos)
xlabel('Eccentricity (mm)')
ylabel('Displaced Location (mm)')
title('New Position Along Radial')
hline =refline(1,0)
hline.Color = 'k'
hline.LineStyle = '--'
axis square
subplot(1,3,3)
plot(diffRadial);
ylabel('Diff (mm)')
title('Diff of Displacement Along Radial')
axis square


end

