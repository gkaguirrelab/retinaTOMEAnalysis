function checkSingleRadial(radialRotDeg,sizeMM,smpPerMM,sectorAngle)


map = makeMapRadial(sizeMM,smpPerMM,sectorAngle);

%Create sample base in MM based on map size and sample rate
sampleBase=0:1/smpPerMM:sizeMM;

rotMap   = imrotate(map,-1.*radialRotDeg,'crop','nearest');

center = round(size(rotMap,1)/2);

radial = rotMap(center,center:end);

dispDiff = diff(radial);

newPos = sampleBase - radial;

diffRadial = diff(newPos);

figure;
subplot(2,2,2)
plot(dispDiff);
ylabel('Diff (mm)')
title('Diff of Displacement Along Radial')
subplot(2,2,1)
plot(sampleBase,radial)
xlabel('Eccentricity (mm)')
ylabel('Displacement (mm)')
title('Displacement Along Radial')
axis square
subplot(2,2,3)
plot(sampleBase,newPos)
xlabel('Eccentricity (mm)')
ylabel('Displaced Location (mm)')
title('New Position Along Radial')
hline =refline(1,0)
hline.Color = 'k'
hline.LineStyle = '--'
axis square
subplot(2,2,4)
plot(diffRadial);
ylabel('Diff (mm)')
title('Diff of New Position Along Radial')
axis square


end

