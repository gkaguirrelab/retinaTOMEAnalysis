
exampleMontage = load(fullfile(dropboxBaseDir,'AOSO_analysis\OCTExplorerExtendedHorizontalData\11028\OD\11028_OD_ExtendedOCT.mat'));
exampleMontageOutlined = imread(fullfile(dropboxBaseDir,'AOSO_analysis\OCTExplorerExtendedHorizontalData\11028\OD\AvgAll_Trans_wManualSmooth_Padded_p15.tif'));

imwrite(exampleMontage.imIn{1},fullfile(saveDir,'fig1','a.png'));
imwrite(exampleMontage.imIn{2},fullfile(saveDir,'fig1','b.png'));
imwrite(exampleMontage.imIn{3},fullfile(saveDir,'fig1','c.png'));
imwrite(exampleMontage.imOut,fullfile(saveDir,'fig1','d.png'));
imwrite(exampleMontageOutlined,fullfile(saveDir,'fig1','f.png'));
