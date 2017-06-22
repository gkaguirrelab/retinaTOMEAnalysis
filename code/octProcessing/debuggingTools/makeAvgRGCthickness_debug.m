% make and average warped rgc thickness map.

%% Set paths and parameters
baseDir  = '~/Dropbox-Aguirre-Brainard-Lab/AOSO_analysis/connectomeRetinaData/';
eye      = 'OD'; % either OD or OS
scans     = {'H','V'};
subjects = {'11015','11018','11028','11043','11050','11051','11052','11053','11055','11056','11057'};
sectorAngle = 3;

%% Load Sample Subject For Sizing 
inFile = fullfile(baseDir,subjects{1},'HeidelbergSpectralisOCT',eye,[subjects{1} '_' eye scans{1}],[subjects{1} '_result.mat']);
%% Generate Displacement Map
dispMap = makeMapFromOCT(inFile,sectorAngle);

randMat = rand(size(dispMap,1),size(dispMap,2),length(subjects)*10)-0.5;
load(inFile)
[rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header);
%% Load and warp all the RGC+ thickness maps
count = 1;
for s = 1:length(subjects)
    for sc = 1:length(scans)
        
        %inFile = fullfile(baseDir,subjects{s},'HeidelbergSpectralisOCT',eye,[subjects{s} '_' eye scans{sc}],[subjects{s} '_result.mat']);
        %load(inFile);
        
        %calculate thickness in mm
%         if strcmp(scans{sc},'H')
%             [rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header);
%         elseif strcmp(scans{sc},'V')
%             %is vertical, rotate by 90 deg to match horizontal
%             [rgcPlus,sampleBaseRadius] = rgcThickness(bd_pts,header);
%              rgcPlus = imrotate(rgcPlus,90);
%         end

        % 
        rgcPlus = randMat(:,:,count);
        
        % generate a sampleBase in mm. The origin (0,0) is at the center
        xAxis=[fliplr(sampleBaseRadius)*-1,sampleBaseRadius(2:end)];
        yAxis=[fliplr(sampleBaseRadius)*-1,sampleBaseRadius(2:end)]';
        sampleBaseX=repmat(xAxis,length(yAxis),1);
        sampleBaseY=repmat(yAxis,1,length(yAxis));
        
        %% Warp Image 
        rgcPlusWarpMat(:,:,count) = warpImage(rgcPlus, dispMap, sampleBaseX, sampleBaseY);
        
%         %% Fill in NaN's
%         rgcPlusWarp(:,:,count) = fillNansInMap(rgcPlusWarp(:,:,count),radMM,smpPerMM,4);
        
        count= count+1;
        clear bd_pts header inFile rgcPlus
    end
end

rgcPlusWarp = nanmean(rgcPlusWarpMat,3);


radMM = max(sampleBaseRadius);
smpPerMM = (length(sampleBaseRadius)-1)./radMM;
diskSize= 3;
warpMap = fillNansInMap(rgcPlusWarp,radMM,smpPerMM,diskSize);


plotResults = 'OFF';

if strcmp(plotResults,'full')
    figure
    mdPt= round(size(warpMap,1)/2);
    subplot(1,2,1)
    plot(-3:1/16:3,warpMap(mdPt,:))
    ylabel('Thickness mm')
    xlabel('Eccentricty mm')
    legend('Vertical Meridian')
    
    subplot(1,2,2)
    plot(-3:1/16:3,warpMap(:,mdPt))
    ylabel('Thickness mm')
    xlabel('Eccentricty mm')
    legend('Horizontal Meridian')
    
    figure
    nasal = test(mdPt,mdPt:end);
    temporal = fliplr(test(mdPt,1:mdPt));
    superior = fliplr(test(1:mdPt,mdPt)');
    inferior = test(mdPt:end,mdPt)';
    figure; hold on
    plot(0:1/16:3,nasal,'g')
    plot(0:1/16:3,temporal,'r')
    plot(0:1/16:3,superior,'b')
    plot(0:1/16:3,inferior,'k')
    legend('inferior','temporal','superior','inferior')
    ylabel('Thickness mm')
    xlabel('Eccentricty mm')
    
    
%     nasal = warpMap(mdPt,mdPt:end);
%     temporal = fliplr(warpMap(mdPt,1:mdPt));
%     superior = fliplr(warpMap(1:mdPt,mdPt)');
%     inferior = warpMap(mdPt:end,mdPt)';
%     ax = figure; hold on
%     ax = plot(1/16:1/16:3,nasal(2:end),'g')
%     ax = plot(1/16:1/16:3,temporal(2:end),'r')
%     ax = plot(1/16:1/16:3,superior(2:end),'b')
%     ax = plot(1/16:1/16:3,inferior(2:end),'k')
%     legend('inferior','temporal','superior','inferior')
%     ylabel('Thickness mm')
%     xlabel('Eccentricty mm')
%     set(ax,'XScale','log');
%     set(ax,'YScale','log');
    
end

