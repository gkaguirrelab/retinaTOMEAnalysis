function fig2d_ExtremeEyes(saveDir)

eye = modelEyeParameters();
a0 = eye.landmarks.fovea.degField;
theta = linspace(0,2*pi,100);
radius=15;
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    rayPath = calcNodalRayFromField(eye,a0+[x y]);
    X(ii,:) = rayPath(:,end);
end
rayPath = calcNodalRayFromField(eye,a0);
plotOpticalSystem(eye,'surfaceAlpha',0.5);
%plot3(X(:,1),X(:,2),X(:,3),'-r');
patch(X(:,1),X(:,2),X(:,3),[1 0 0])
axis off
print('~/Desktop/Eye30DegPatch.png','-dpng','-r600')

bigEyeRadiusMm = mean(vecnorm((X-mean(X))')');



eye = modelEyeParameters('axialLength',21.79,'sphericalAmetropia',3.875);
a0 = eye.landmarks.fovea.degField;
theta = linspace(0,2*pi,20);
radius=15;
for ii=1:length(theta)
    x=radius*cos(theta(ii));
    y=radius*sin(theta(ii));
    rayPath = calcNodalRayFromField(eye,a0+[x y]);
    X(ii,:) = rayPath(:,end);
end
rayPath = calcNodalRayFromField(eye,a0);
plotOpticalSystem(eye,'newFigure',false);
plot3(X(:,1),X(:,2),X(:,3),'.r');

smallEyeRadiusMm = mean(vecnorm((X-mean(X))')');

xlim([-30 10]);
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