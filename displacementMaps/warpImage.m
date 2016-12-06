function [ outputImage ] = warpImage( sourceImage, displacementMap, sampleBaseX, sampleBaseY )

% Each point in the displacementMap is the magniude of a vector that is
% directed towards the origin of the image. For each point in the
% sourceImage, we calculate the distance (radius) of the point from the
% origin. The new distance (radius) for that point will be the current
% distance minus the magnitude of displacement. The new distance is then
% decomposed into the X and Y axis components, using the property that the
% ratio of the X and Y distances will be the same before and after
% displacement.


% Calculate the radius for the displaced point
radiusImage=sqrt(sampleBaseX.^2 + sampleBaseY.^2);
newRadiusImage=radiusImage-displacementMap;


% Caculate the ratioImage. To handle infinite values at X=0,
%  replace these with realmax
ratioImage=sampleBaseY./sampleBaseX;
infIdx = find(isinf(ratioImage));
if ~isempty(infIdx)
    ratioImage(infIdx)=realmax;
end

% Decompose the radius of the displaced point into the X and Y position
newX=sqrt(newRadiusImage.^2 ./ (ratioImage.^2 + 1));
newY=newX.*ratioImage;

% Handle the special case of the vertical meridian for newY
vertMeridianIdx=find(sampleBaseX==0);
if ~isempty(radiusImage)
    newY(vertMeridianIdx)=newRadiusImage(vertMeridianIdx);
    newY(vertMeridianIdx)=newY(vertMeridianIdx) .* sign(sampleBaseY(vertMeridianIdx));
end

% Handle the special case of the origin for newX
nanIdx=find(isnan(newX));
if ~isempty(nanIdx)
    newX(nanIdx)=0;
end

% Correct sign for coordinate positions in X
newX=newX.* sign(sampleBaseX);

% Correct sign for coordinate positions in Y
negativeXidx=find(sampleBaseX<0);
if ~isempty(negativeXidx)
    newY(negativeXidx)=newY(negativeXidx).* -1;
end

% restore nans to the newX and newY where displacement is undefined
nanIdx=find(isnan(displacementMap));
if ~isempty(nanIdx)
    newX(nanIdx)=nan;
    newY(nanIdx)=nan;
end

% Create an output image space with double the resolution of the sampleBase
hiresBaseX=imresize(sampleBaseX,[size(sampleBaseX,1)*2 size(sampleBaseX,2)*2],'bilinear');
hiresBaseY=imresize(sampleBaseY,[size(sampleBaseY,1)*2 size(sampleBaseY,2)*2],'bilinear');

% Run through each point in the sourceImage and add it to the
% hiresOutputImage at the closest position identified by newX, newY
xCoordRow=hiresBaseX(1,:);
yCoordRow=hiresBaseY(:,1);
hiresOutputImage=hiresBaseX.*0;

for ii=1:size(sourceImage,1)
    for jj=1:size(sourceImage,2)
        if ~isnan(newX(ii,jj))
            diffX=abs(xCoordRow-newX(ii,jj));
            id_x=find(diffX==min(diffX));
            id_x=id_x(1);
            
            diffY=abs(yCoordRow-newY(ii,jj));
            id_y=find(diffY==min(diffY));
            id_y=id_y(1);
            hiresOutputImage(id_x,id_y)=hiresOutputImage(id_x,id_y)+sourceImage(ii,jj);
        end
    end
end

% Wherever the output image has a value of zero, change it to nan
zeroIdx=find(hiresOutputImage==0);
if ~isempty(zeroIdx)
    hiresOutputImage(zeroIdx)=nan;
end

% Downsample the outputImage to the original size by local averaging,
% avoiding the nans
fac = 0.5;
dx = 1./fac;
[r,c] = ndgrid(1:size(hiresOutputImage,1), 1:size(hiresOutputImage,2));
[n, ibin] = histc(r(:), 0.5:dx:size(hiresOutputImage,1)+0.5);
[n, jbin] = histc(c(:), 0.5:dx:size(hiresOutputImage,2)+0.5);
nr = max(ibin);
nc = max(jbin);
idx = sub2ind([nr nc], ibin, jbin);
outputImage = accumarray(idx, hiresOutputImage(:), [nr*nc 1], @nanmean);
outputImage = reshape(outputImage, nr, nc);

% % Plot the results
% figure
% imagesc(sourceImage)
% figure
% imagesc(outputImage)

% Create x,y,z vectors of the outputImage to allow spline fitting
% valIdx=find(~isnan(outputImage));
% x=sampleBaseX(valIdx);
% y=sampleBaseY(valIdx);
% z=outputImage(valIdx);
% sf=fit([x, y], z,'thinplateinterp');
% figure
% plot(sf,[x,y],z);

end % warp function

