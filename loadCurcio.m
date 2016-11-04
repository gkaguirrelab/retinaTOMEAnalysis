function [dRGC] = loadCurcio(fname,radDeg,smpPerDeg)

fileID = fopen(fname,'r');
C = textscan(fileID,'%f %s %f %f %f');
fclose(fileID);

R = 11.459;

%[X Y] = pol2cart(deg2rad(C{4}),deg2rad(C{3}));
[X,Y] = pol2cart(C{4},C{3});
F = scatteredInterpolant([X Y],C{5});
F.Method = 'linear';
F.ExtrapolationMethod='linear';


%lin = -round(14):stepSize:round(14);
lin = linspace(-radDeg,radDeg,2.*radDeg*smpPerDeg);

[x,y] = meshgrid(lin,lin);
c=((x.^2+y.^2).^.5<=radDeg);

x(c~=1) = nan;
y(c~=1) = nan;
dRGC = F(x,y);


end

