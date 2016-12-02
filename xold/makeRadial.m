function [pol,ecc] = makeRadial(X,Y)

XnumPoints = (2*X);
YnumPoints = (2*Y);

Xpoints = linspace(-1*(X),X,XnumPoints);
Ypoints = linspace(-1*(Y),Y,YnumPoints);

[X1,Y1] = meshgrid(Xpoints,Ypoints);

rad = X1.^2+Y1.^2;

[pol,ecc] = cart2pol(X1,Y1);

pol = flipud(imrotate(rad2deg(pol)+180,180));

end