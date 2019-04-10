inDirVolume = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\volumeMapsBySubject';
inDirThickness = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\averageThicknessMapsBySubject';
outDir = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\thicknessVsVolumeComparison';
subIDs = dir(fullfile(inDirVolume,'1*'));

%We're going to look at some measurements, this describes
measurements = zeros(length(subIDs),9);%this will be out ouput matrix
%This describes each column in the measurements
header = {'subject ID', 'RGCIPL mean thickness', 'RGCIPL mean volume', ...
    'RNFL mean thickness', 'RNFL mean volume', ...
    'OPL mean thickness', 'OPL mean volume', ...
    'Total Retina mean thickness', 'Total Retina mean volume'};


for L = 1:4 %L controls which layer we're looking
    overlap = [];
    
    %first we find the overlaping region across all subjects
    for i = 1:length(subIDs)
        LoadthicknessMap=load(fullfile(inDirThickness,subIDs(i).name,[subIDs(i).name '_averageMaps.mat']));
        
        switch L
            case 1
                thicknessMap = LoadthicknessMap.averageMaps.RGCIPL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_RGCIPL_volumeMaps.mat']);
                volumeMaps = load(loadname);
                
            case 2
                thicknessMap = LoadthicknessMap.averageMaps.RNFL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_RNFL_volumeMaps.mat']);
                volumeMaps = load(loadname);
            case 3
                thicknessMap = LoadthicknessMap.averageMaps.OPL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_OPL_volumeMaps.mat']);
                volumeMaps = load(loadname);
            case 4
                thicknessMap = LoadthicknessMap.averageMaps.TotalRetina;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_TotalRetina_volumeMaps.mat']);
                volumeMaps = load(loadname);
        end
        
        if(isempty(overlap))
            overlap = ones(size(volumeMaps.volumeMap_mmCubed));
        end
        
        overlap= overlap & ~isnan(volumeMaps.volumeMap_mmCubed) &  ~isnan(thicknessMap);
    end
    
    %writeout the overlap
    
    switch L
        case 1
            savename = fullfile(outDir, 'RGCIPL_overlapMap.mat');
        case 2
            savename = fullfile(outDir, 'RNFL_overlapMap.mat');
        case 3
            savename = fullfile(outDir, 'OPL_overlapMap.mat');
        case 4
            savename = fullfile(outDir, 'TotalRetina_overlapMap.mat');
    end
    save(savename,'overlap');
    
    
    %now that we've got the overlap, we go back through and calculate the
    %mean for each subject across the overlap
    
    for i = 1:length(subIDs)
        LoadthicknessMap=load(fullfile(inDirThickness,subIDs(i).name,[subIDs(i).name '_averageMaps.mat']));
        
        switch L
            case 1
                thicknessMap = LoadthicknessMap.averageMaps.RGCIPL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_RGCIPL_volumeMaps.mat']);
                volumeMaps = load(loadname);
                
            case 2
                thicknessMap = LoadthicknessMap.averageMaps.RNFL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_RNFL_volumeMaps.mat']);
                volumeMaps = load(loadname);
            case 3
                thicknessMap = LoadthicknessMap.averageMaps.OPL;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_OPL_volumeMaps.mat']);
                volumeMaps = load(loadname);
            case 4
                thicknessMap = LoadthicknessMap.averageMaps.TotalRetina;
                loadname = fullfile(inDirVolume, subIDs(i).name, [subIDs(i).name '_TotalRetina_volumeMaps.mat']);
                volumeMaps = load(loadname);
        end
        
        measurements(i,1) = str2double(subIDs(i).name);
        measurements(i,2*L) = mean(thicknessMap(overlap));
        measurements(i,2*L+1) =  mean(volumeMaps.volumeMap_mmCubed(overlap));
    end
    
end


filename = fullfile(outDir,'meanThicknessAndVolumes.xlsx');
xlswrite(filename,header,1,'A1');
xlswrite(filename,measurements,1,'A2');