function AnalyzeGCIPExtendedOCT(dataDir)
%Purpose does some preliminary data analysis of the extended OCT
%segmentations
%find all the manual boundaries in the directories
allFiles = subdir(fullfile(dataDir,'*ManualBoundaries.mat'));
N =  length(allFiles);

%set figures
figure(1)
clf
set (gca,'YDir','reverse')
hold on

figure(2)
clf
hold on

figure(3)
clf
hold on

%go through all the subjects
for i = 1:N

    if(i == 42)%We're skipping some subjects because of errors in the segmentation
        continue
    end
    
    %load the boundaries
    load(allFiles(i).name);
    randcolor = [rand(1) rand(1) rand(1)];%set a random color for this subject

    %if OS, then we flip
    if(~isempty(strfind(allFiles(i).name,'\OS\')))
       flipSign = -1;
       segLeft = 2;
       segRight = 3;

    else
       flipSign = 1; 
       segLeft = 1;
       segRight = 2;
    end
    
    
    %find the fovea location
       [leftLoc_x, leftLoc_ind] = max(boundariesSmooth(segLeft).GC_IP(:,1));
       leftLoc_y = boundariesSmooth(segLeft).GC_IP(leftLoc_ind,2);
       [rightLoc_x, rightLoc_ind] = min(boundariesSmooth(segRight).GC_IP(:,1));
       rightLoc_y = boundariesSmooth(segRight).GC_IP(rightLoc_ind,2);
       foveaLoc_x = flipSign*mean([leftLoc_x rightLoc_x]);
       foveaLoc_y = mean([leftLoc_y rightLoc_y]);

    %Plot 1 show all the boundaries we found across all subjects
    figure(1)
    for c = 1:3
    plot(flipSign*boundariesSmooth(c).GC_Inner(:,1)-foveaLoc_x,boundariesSmooth(c).GC_Inner(:,2)-foveaLoc_y,'color',randcolor(:))
    plot(flipSign*boundariesSmooth(c).GC_IP(:,1)-foveaLoc_x,boundariesSmooth(c).GC_IP(:,2)-foveaLoc_y,'color',randcolor(:))
    plot(flipSign*boundariesSmooth(c).IP_Outer(:,1)-foveaLoc_x,boundariesSmooth(c).IP_Outer(:,2)-foveaLoc_y,'color',randcolor(:))
    end

    %Plot 2 calculate the minimum distance from GCIP boundary to inner GC
    %boundary
    figure(2)
    for c = 1:3
        plot(flipSign*boundariesSmooth(c).GC_IP(:,1)-foveaLoc_x,boundaryDist(boundariesSmooth(c).GC_IP,boundariesSmooth(c).GC_Inner,1),'color',randcolor(:))
    end
    
    %Plot 3 calculate the minimum distance from GCIP boundary to outer GC
    %boundary

    figure(3)
    for c = 1:3
        plot(flipSign*boundariesSmooth(c).GC_IP(:,1)-foveaLoc_x,boundaryDist(boundariesSmooth(c).GC_IP,boundariesSmooth(c).IP_Outer,1),'color',randcolor(:))
    end
      
end
