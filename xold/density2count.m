function [countDf,countRGC] = density2count(Df,RGC,areaVec)

midPt = round(size(Df,1)/2);


dfVals = Df(midPt,midPt:end);
rgcVals = RGC(midPt,midPt:end);

countDf=dfVals.*areaVec;
countRGC=rgcVals.*areaVec;
end