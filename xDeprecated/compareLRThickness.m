dataDir='C:\Users\dontm\Dropbox (Aguirre-Brainard Lab)\AOSO_analysis\OCTExplorerExtendedHorizontalData';
allFiles = dir(fullfile(dataDir,'1*'));

N = length(allFiles)
figure(1)
clf
hold on
figure(2)
clf
hold on
for n = 1:N
    figure(1)
    clf
    hold on
    figure(2)
    clf
    hold on
    figure(3)
    hold on

    OD_FILE = fullfile(dataDir,allFiles(n).name,'OD','ManualBoundaries.mat');
    OS_FILE = fullfile(dataDir,allFiles(n).name,'OS','ManualBoundaries.mat');
    if(~isfile(OD_FILE) || ~isfile(OS_FILE))
        continue
    end
    
    randcolor = [rand(1) rand(1) rand(1)];%set a random color for this subject
    
    for c = 1:3
        %OD
        cc = c;
        load(OD_FILE);
        flipSign = 1;
        segLeft = 1;
        segRight = 2;
        
        %find the fovea location
        [leftLoc_x, leftLoc_ind] = max(boundariesSmooth(segLeft).GC_IP(:,1));
        leftLoc_y = boundariesSmooth(segLeft).GC_IP(leftLoc_ind,2);
        [rightLoc_x, rightLoc_ind] = min(boundariesSmooth(segRight).GC_IP(:,1));
        rightLoc_y = boundariesSmooth(segRight).GC_IP(rightLoc_ind,2);
        foveaLoc_x = flipSign*mean([leftLoc_x rightLoc_x]);
        foveaLoc_y = mean([leftLoc_y rightLoc_y]);
        OD_xLoc = flipSign*boundariesSmooth(cc).GC_IP(:,1)-foveaLoc_x;
        OD_GCthickness = boundaryDist(boundariesSmooth(cc).GC_IP,boundariesSmooth(cc).GC_Inner,1);
        OD_IPthickness = boundaryDist(boundariesSmooth(cc).GC_IP,boundariesSmooth(cc).IP_Outer,1);
        figure(1)
        plot(OD_xLoc,OD_IPthickness,'color',randcolor(:))
        
        %OS
        cc = 4-c;
        load(OS_FILE);
        flipSign = -1;
        segLeft = 2;
        segRight = 3;
        %find the fovea location
        [leftLoc_x, leftLoc_ind] = max(boundariesSmooth(segLeft).GC_IP(:,1));
        leftLoc_y = boundariesSmooth(segLeft).GC_IP(leftLoc_ind,2);
        [rightLoc_x, rightLoc_ind] = min(boundariesSmooth(segRight).GC_IP(:,1));
        rightLoc_y = boundariesSmooth(segRight).GC_IP(rightLoc_ind,2);
        foveaLoc_x = flipSign*mean([leftLoc_x rightLoc_x]);
        foveaLoc_y = mean([leftLoc_y rightLoc_y]);
        OS_xLoc = flipSign*boundariesSmooth(cc).GC_IP(:,1)-foveaLoc_x;
        OS_GCthickness = boundaryDist(boundariesSmooth(cc).GC_IP,boundariesSmooth(cc).GC_Inner,1);
        OS_IPthickness = boundaryDist(boundariesSmooth(cc).GC_IP,boundariesSmooth(cc).IP_Outer,1);
        figure(1)
        
        plot(OS_xLoc,OS_IPthickness,'color',randcolor(:))
        eyeDiff = boundaryDist([OS_xLoc OS_IPthickness], [OD_xLoc OD_IPthickness],0);
        
        figure(2)
        plot(OS_xLoc,eyeDiff,'color',randcolor(:))

        figure(3)
        plot(OS_xLoc,eyeDiff,'color',randcolor(:))
    end
    

    pause
end

