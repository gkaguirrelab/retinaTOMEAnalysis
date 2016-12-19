function diffRadial = checkDispMatrix(radialDegree,sizeMM,smpPerMM,sectorAngle)


map = makeMap(sizeMM,smpPerMM,sectorAngle);


%Create sample base in MM based on map size and sample rate
sampleBase=0:1/smpPerMM:sizeMM;
sampleBaseSurf = repmat(sampleBase,size(map,1),1);

newPos = sampleBaseSurf - map;

diffRadial = diff(newPos,1,2);


rowNum = find(min(abs((0:sectorAngle:360)-radialDegree)))

% figure;
% surf(diffRadial);

figure;
subplot(2,2,1)
imagesc(map)
subplot(2,2,2)
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
title('Diff of Displacement Along Radial')
axis square


end

