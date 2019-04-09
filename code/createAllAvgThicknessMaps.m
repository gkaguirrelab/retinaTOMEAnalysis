load('C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\averageThicknessMapsBySubject\everyThicknessMap.mat');
outdir = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\averageThicknessMapsBySubject';

for i = 1:length(subjectList)
    i
    subdir = fullfile(outdir,subjectList{i});
    mkdir(subdir);
    
    %load up everyThicknessMap(layernames,subjectInd,OS/OD,XN,YN)
    %layernames:1 - RGC+IPL',2-'RNFL',3-'OPL',4-'total
    %OS/OD - 1-OS, 2-OD
    
    %save the average of the OS (which has already been mirroed) and OD
    averageMaps.RGCIPL = squeeze(mean(everyThicknessMap(1,i,:,:,:),3));
    averageMaps.RNFL = squeeze(mean(everyThicknessMap(2,i,:,:,:),3));
    averageMaps.OPL = squeeze(mean(everyThicknessMap(3,i,:,:,:),3));
    averageMaps.TotalRetina = squeeze(mean(everyThicknessMap(4,i,:,:,:),3));

    savename = fullfile(subdir,[subjectList{i} '_averageMaps.mat']);
    save(savename,'averageMaps');
   
end

    