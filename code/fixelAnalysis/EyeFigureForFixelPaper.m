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