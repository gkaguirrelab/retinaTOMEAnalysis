function dispMap = fillNansInMap(nanDispMap,radMM,smpPerMM,diskSize)

 [id1,id2]= find(isnan(nanDispMap));
 dispMap = nanDispMap;
 for i = 1:length(id1)
     dot = zeros(size(nanDispMap));
     dot(id1(i),id2(i)) = 1;
     se = strel('disk',diskSize,0);
     mask = imdilate(dot,se);
     dispMap(id1(i),id2(i)) = nanmean(nanDispMap(logical(mask)));
 end
     
 %% mask the output
 line = -radMM:1/smpPerMM:radMM;
 [X1,X2] = meshgrid(line,line);
 mask = sqrt(X1.^2 + X2.^2) <= radMM;
 
 dispMap(~mask) = nan;
end
 