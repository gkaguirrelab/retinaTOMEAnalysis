% make and average warped rgc thickness map.

baseDir  = '~/Dropbox-Aguirre-Brainard-Lab/AOSO_analysis/connectomeRetinaData/'; 
eye      = 'OD'; % either OD or OS
scans     = {'H','V'};
subjects = {'11015','11018','11028','11043','11050','11051','11052','11053','11055','11056','11057'};

for s = 1:length(subjects)
    for sc = 1:length(scans)
    
        inFile = fullfile(baseDir,subjects{s},'HeidelbergSpectralisOCT',eye,[subjects{s} '_' eye scans{sc}],[subjects{s} '_result.mat']);
        load(inFile);
        if strcmp(scans{sc},'H')
            [rgcPlus(:,:,s+sc-1),sampleBaseRadius] = rgcThickness(bd_pts,header);
        elseif strcmp(scans{sc},'V')
            [rgcPlus(:,:,s+sc-1),sampleBaseRadius] = rgcThickness(bd_pts,header);
             rgcPlus(:,:,s+sc-1) = imrotate(rgcPlus(:,:,s+sc-1),90);
        end
    end
end

        
