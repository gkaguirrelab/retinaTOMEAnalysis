radDeg = 20;
smpPerDeg=2;
sectorAngle = 6;

%% Generate the Rential Ganglion Cell Density
RGCdensity= densityRGC(radDeg,smpPerDeg,'OFF');

%% Generate the Recptive Field Density
RFdensity = densityRf(radDeg,smpPerDeg,'OFF');

areaPerSeg= calcSegSize(radDeg,smpPerDeg,sectorAngle);

[countRF,countRGC] = density2count(RFdensity,RGCdensity,areaPerSeg);

countRFsum = cumsum(countRF);
countRGCsum = cumsum(countRGC);

validatePlots(RGCdensity,RFdensity,countRF,countRGC,radDeg,smpPerDeg);
