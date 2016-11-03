
radDeg = 14;
smpPerDeg = 4;
fname = '/Users/michael/Documents/MATLAB/octAnalysisForTOME/curcio_data.txt';


[line,pol,ecc] = createGrid(radDeg,smpPerDeg);

Df = densityRf(ecc,pol);

[dRGC] = loadCurcio(fname,radDeg,smpPerDeg);


figure
subplot(2,1,1);
imagesc(Df)
title('RGC receptive feild Density')
axis square
subplot(2,1,2);
imagesc(dRGC)
title('RGC Density')
axis square


