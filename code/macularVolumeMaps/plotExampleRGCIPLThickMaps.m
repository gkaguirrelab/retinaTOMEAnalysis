dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
dataDir = fullfile(dropboxBaseDir,'AOSO_analysis','averageThicknessMapsBySubject');

subjectList = {'11050','11083','11100'};
figA = figure();

for ii=1:length(subjectList)
    fileName = fullfile(dataDir,subjectList{ii},[subjectList{ii} '_averageMaps.mat']);
    load(fileName,'averageMaps');
    thisMap = averageMaps.RGCIPL;
    imageSize = size(thisMap);
    % Convert from microns to mm
    thisMap = thisMap./1000;
    
    zMax = 0.12;
    supportDeg = linspace(-15,15,imageSize(1));
    [xMesh,yMesh] = meshgrid(supportDeg,supportDeg);
    
    
    figure(figA);
    subplot(2,3,ii)
    mesh(xMesh,yMesh,thisMap);
    caxis([0 zMax]);
    xlim([-15 15]);
    ylim([-15 15]);
    zlim([0 zMax]);
    zlabel('tissue thickness [mm]');
    title(['Subject ' subjectList{ii}] )
    axis square
end
figSaveDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/pca';
filename = fullfile(figSaveDir,'rgciplThickMapExamples');
vecrast(figA, filename, 600, 'top', 'pdf')
close(figA)