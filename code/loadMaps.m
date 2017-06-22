function [rgcPlus,dispMap, sampleBaseX, sampleBaseY] = loadMaps(octSegFile,sectorAngle)
% load the output of Arua Tools and get the dimensions generates a matching
% sized dispalcemnt map.

%% preset sample data if no data input defined
if isempty(octSegFile)
    base_dir = path2oct;
    octSegFile = [base_dir '/octProcessing/sample_data/11015_OD_Horizontal/11015_1_21682_result.mat'];
end
if isempty(sectorAngle) 
    sectorAngle=6;
end

% load OTC segmentation
load(octSegFile)
% reate RGC+ thickness and sampleBase
[rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header);
% Turn sambleBase into intputs for makeMap
radMM = max(sampleBaseRadius);
smpPerMM = (length(sampleBaseRadius)-1)./radMM;
% generate displacement map
dispMap = makeMap(radMM,smpPerMM,sectorAngle);
% generate a sampleBase in mm. The origin (0,0) is at the center
xAxis=[fliplr(sampleBaseRadius)*-1,sampleBaseRadius(2:end)];
yAxis=[fliplr(sampleBaseRadius)*-1,sampleBaseRadius(2:end)]';
sampleBaseX=repmat(xAxis,length(yAxis),1);
sampleBaseY=repmat(yAxis,1,length(yAxis));

end

