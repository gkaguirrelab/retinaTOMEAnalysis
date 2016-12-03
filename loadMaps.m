function [rgcPlus,dispMap] = loadMaps(octSegFile,sectorAngle)
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
% vreate RGC+ thickness and sampleBase
[rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header);
% Turn sambleBase into intputs for makeMap
radMM = max(sampleBaseRadius);
smpPerMM = (length(sampleBaseRadius)-1)./radMM;
% generate displacement map
dispMap = makeMap(radMM,smpPerMM,sectorAngle);

end

