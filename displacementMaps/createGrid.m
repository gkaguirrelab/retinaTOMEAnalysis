function [line,pol,ecc] = createGrid(radMM,smpPerMM)

line = linspace(-radMM,radMM,2.*radMM*smpPerMM+1);


Xpoints = line;
Ypoints = line;

[X1,Y1] = meshgrid(Xpoints,Ypoints);


[pol,ecc] = cart2pol(X1,Y1);

pol = flipud(imrotate(rad2deg(pol)+180,180));

end