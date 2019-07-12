function Dist = boundaryDist(B1, B2, DistType)
%Calculates the minimum distance between each point in Curve B1 and Curve B2
%Inputs
%B1 and B2 - Input boundaries
%DistType: 0=vertical Distance, 1=Minimum Distance
%Outputs
%Dist - returns the distance at each location in B1 to B2

if (DistType == 0)
    [Idx,Dist] = knnsearch(Surf2(:,1), Surf1(:,1));%find closes vertical x location
    Dist = abs(B1(:,2) - B2(Idx,2));
elseif (DistType == 1)
    [Idx,Dist] = knnsearch(B2, B1);%we just use NN in this case
end
