function displace = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg)

%%% get the displacement map working 

displaceMap = nan(size(retData,1),size(retData,2),upper(360/sectorAngle));
for i = 1:upper(360/sectorAngle) 
    [rotDf,rotRGC] = rotAndCrop(retData(:,:,4),retData(:,:,1),rotDeg);
    areaPerSeg= calcSegSize(radDeg,smpPerDeg,sectorAngle);
    degVec = 0:1/smpPerDeg:radDeg;
    [countDf,countRGC] = density2count(rotDf,rotRGC,areaPerSeg);
    [DfFit] = fit(degVec',cumsum(countDf)','smoothingspline','Exclude', find(isnan(countDf)),'SmoothingParam', 1);
    [RGCFit] = fit(degVec',cumsum(countRGC)','smoothingspline','Exclude', find(isnan(countDf)),'SmoothingParam', 1);



    rgcVal = RGCFit(degVec);
    eccDf = interp1(DfFit(degVec),degVec,rgcVal,'spline');
    displace = abs(degVec'-eccDf);

    %displace(find(displace == min(displace(2:end)))+1:end) = nan;
    
    %displaceMap = displaceMap(round(size
end


end
