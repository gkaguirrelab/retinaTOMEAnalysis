%load data
OCTIn = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\thicknessVsVolumeComparison\meanThicknessAndVolumes.csv';
MRIIn = 'C:\Users\dontm\Dropbox (Personal)\Research\Projects\retinaTOMEAnalysis\data\visualPathwayAnatMeasures.csv';
outFile = 'C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\thicknessVsVolumeComparison\OCTvsMR.xls';

OCTData =  csvread(OCTIn);

[num,txt,raw] = xlsread(MRIIn);
MRIData = num;

[Lia, LocB] = ismember(MRIData(:,1),OCTData(:,1),'rows');



MRIDataToOCT = zeros(size(OCTData,1),size(MRIData,2));

MRIDataToOCT(LocB,:) = MRIData;

xlswrite(outFile,OCTData,1,'A1')

xlswrite(outFile,MRIDataToOCT,1,'P1');
