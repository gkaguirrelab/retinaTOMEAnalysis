function displace = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg)

%%% get the displacement map working

[rotDf,rotRGC] = rotAndCrop(retData(:,:,4),retData(:,:,1),rotDeg);
areaPerSeg= calcSegSize(radDeg,smpPerDeg,sectorAngle);
degVec = 1/smpPerDeg:1/smpPerDeg:radDeg+1/smpPerDeg;
[countDf,countRGC] = density2count(rotDf,rotRGC,areaPerSeg);
[DfFit] = fit(degVec',cumsum(countDf)','smoothingspline','Exclude', find(isnan(countDf)),'SmoothingParam', 1);
[RGCFit] = fit(degVec',cumsum(countRGC)','smoothingspline','Exclude', find(isnan(countDf)),'SmoothingParam', 1);



rgcVal = RGCFit(degVec);
eccDf = interp1(DfFit(degVec),degVec,rgcVal,'spline');
displace = abs(degVec'-eccDf);




end 
