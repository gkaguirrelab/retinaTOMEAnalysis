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

bigEyeRadiusMm = mean(vecnorm((X-mean(X))')');


sceneGeometry = createSceneGeometry('axialLength',21.79,'sphericalAmetropia',3.875);
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    [~,X(ii,:)] = calcRetinaFieldPoint( sceneGeometry.eye, a0+[x y 0]);
end
plotOpticalSystem('newFigure',false,'surfaceSet',sceneGeometry.refraction.retinaToStop);
plot3(X(:,1),X(:,2),X(:,3),'.r');

smallEyeRadiusMm = mean(vecnorm((X-mean(X))')');

title('TOME_3031, 3043, SR=-10, 4, axLength=27.57, 21.79');
set(gca,'color','none')
drawnow
filename =fullfile(saveDir,'fig2','d_eye');
vecrast(gcf, filename, 600, 'bottom', 'pdf')
close(gcf)


fig = figure;
fig.Renderer='Painters';
thickA = 0.5;

[X,Y,Z] = cylinder(bigEyeRadiusMm,100);
Z = Z*thickA;

h = surf(X,Y,Z,'facecolor','b','LineStyle','none');
h.FaceAlpha = 0.25;
hold on
axis equal
h = fill3(X(1,:),Y(1,:),Z(1,:),'-b');
h.FaceAlpha = 0.25;
h = fill3(X(2,:),Y(2,:),Z(2,:),'-b');
h.FaceAlpha = 0.25;

V = pi * bigEyeRadiusMm^2 * thickA;
thickB = V / (pi * smallEyeRadiusMm^2);

[X,Y,Z] = cylinder(smallEyeRadiusMm,100);
Z = Z*thickB;

h = surf(X,Y,Z,'facecolor','r','LineStyle','none');
h.FaceAlpha = 0.25;
hold on
axis equal
%h = fill3(X(1,:),Y(1,:),Z(1,:),'r');
%h.FaceAlpha = 0.25;
h = fill3(X(2,:),Y(2,:),Z(2,:),'-r');
h.FaceAlpha = 0.25;

view(-36,26)
axis off
title('retinal cyliners');
set(gca,'color','none')
drawnow
filename =fullfile(saveDir,'fig2','d_cylinder.pdf');

setTightFig
saveas(fig,filename);



end