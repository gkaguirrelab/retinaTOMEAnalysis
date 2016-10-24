function radPoints = makeRadial(X,Y)


XnumPoints = (2*X);
YnumPoints = (2*Y);

Xpoints = linspace(-1*(X),X,XnumPoints);
Ypoints = linspace(-1*(Y),Y,YnumPoints);

[X1,Y1] = meshgrid(Xpoints,Ypoints);

rad = X1.^2+Y1.^2;

%% scale by max value
  
if Y<=X
    radPoints = max([X Y]).*(rad./max(rad(:,Y)));
else
    radPoints = max([X Y]).*(rad./max(rad(X,:)));
end

end