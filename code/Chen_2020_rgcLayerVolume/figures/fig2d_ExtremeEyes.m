function fig2d_ExtremeEyes(saveDir)

sceneGeometry = createSceneGeometry('axialLength',27.57,'sphericalAmetropia',-10.25);
a0 = [5.45 2.5 0];
theta = linspace(0,2*pi,20);
radius=15;
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    [~,X(ii,:)] = calcRetinaFieldPoint( sceneGeometry.eye, a0+[x y 0]);
end
G = calcRetinaFieldPoint( sceneGeometry.eye, a0);
[outputRay,rayPath] = calcNodalRay(sceneGeometry.eye,G);
plotOpticalSystem('surfaceSet',sceneGeometry.refraction.retinaToCamera,'addLighting',true,'rayPath',rayPath,'outputRay',outputRay);
plot3(X(:,1),X(:,2),X(:,3),'.b');

sceneGeometry = createSceneGeometry('axialLength',21.79,'sphericalAmetropia',3.875);
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    [~,X(ii,:)] = calcRetinaFieldPoint( sceneGeometry.eye, a0+[x y 0]);
end
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.retinaToStop);
plot3(X(:,1),X(:,2),X(:,3),'.r');


title('TOME_3031, 3043, SR=-10, 4, axLength=27.57, 21.79');
set(gca,'color','none')
drawnow
filename =fullfile(saveDir,'fig2','d.pdf');
vecrast(gcf, filename, 600, 'bottom', 'pdf')
close(gcf)

end