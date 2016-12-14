function dispMap = fillNansInMap(nanDispMap,radMM,smpPerMM)

[X,Y] = meshgrid(1:size(nanDispMap,1),1:size(nanDispMap,2));

%// identify indices valid for the 3 matrix 
idxgood=~(isnan(X) | isnan(Y) | isnan(nanDispMap)); 

%// re-interpolate scattered data (only valid indices) over the "uniform" grid
dispMap = griddata( X(idxgood),Y(idxgood),nanDispMap(idxgood), X, Y ,'cubic') ;

%  [id1,id2]= find(nanDispMap == nan);
%  dispMap = nanDispMap;
%  for i = 1:length(id1)
%      dot = zeros(size(nanDispMap));
%      dot(id1(i),id2(i)) = 1;
%      se = strel('disk',diskSize,0);
%      mask = imdilate(dot,se);
%      dispMap(id1(i),id2(i)) = nanmean(nanDispMap(logical(mask)));
%  end
     
 %% mask the output
 line = -radMM:1/smpPerMM:radMM;
 [X1,X2] = meshgrid(line,line);
 mask = sqrt(X1.^2 + X2.^2) <= radMM;
 
 dispMap(~mask) = nan;
end
 