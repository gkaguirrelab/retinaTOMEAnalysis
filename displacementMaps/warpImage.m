function [ outputImage ] = warpImage( sourceImage, displacementMap, sampleBaseX, sampleBaseY )

newY=sampleBaseY*0;
ratioImage=sampleBaseY./sampleBaseX;
radiusImage=sqrt(sampleBaseX.^2 + sampleBaseY.^2);
newRadiusImage=radiusImage-displacementMap;
newX=sqrt(newRadiusImage.^2 ./ (ratioImage.^2 + 1));
newY=newX.*ratioImage;


end

