
% This is the mmPerDeg at the ellipsoidal pole of the vitreous chamber
mmPerDeg = @(axialLength) (0.0165.*axialLength)-0.1070;


sceneGeometry = createSceneGeometry('axialLength',27.57,'sphericalAmetropia',-10.25);
a0 = [5.45 2.5 0];
theta = linspace(0,2*pi,20);
radius=15;%just an example
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    [~,X(ii,:)] = findRetinaFieldPoint( sceneGeometry.eye, a0+[x y 0]);
end
G = findRetinaFieldPoint( sceneGeometry.eye, a0);
[outputRay,rayPath] = calcNodalRay(sceneGeometry.eye,G);
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
plot3(X(:,1),X(:,2),X(:,3),'.r');

sceneGeometry = createSceneGeometry('axialLength',21.79,'sphericalAmetropia',3.875);
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    [~,X(ii,:)] = findRetinaFieldPoint( sceneGeometry.eye, a0+[x y 0]);
end
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.retinaToLens);
plot3(X(:,1),X(:,2),X(:,3),'.r');



title(['TOME_3031, 3043, SR=-10, 4, axLength=27.57, 21.79, mm/deg=' num2str(mmPerDeg(27.57)) ', ' num2str(mmPerDeg(21.79))]);
set(gca,'color','none')
drawnow
filename = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/modelEye';
vecrast(gcf, filename, 600, 'bottom', 'pdf')
close(gcf)
