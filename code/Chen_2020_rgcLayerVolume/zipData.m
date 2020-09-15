inDir='D:\Min\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\OCTSingleVerticalData'
outDir='C:\Users\dontm\Downloads\temp\New folder (2)'

allFiles = dir(inDir);

allFiles(1)=[];
allFiles(1)=[];

for i = 1:length(allFiles)
   
    inname=fullfile(inDir,allFiles(i).name)
    outname=fullfile(outDir,[allFiles(i).name '.zip'])
    

    zip(outname,inname)
    
    
end