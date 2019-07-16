function smoothManualSegGCIPExtOCT(dataDir)
%Purpose: Find all manual segmentations in dataDir with Horizontal Extended
%OCT data and calculate and smooth the the manual segmentation boundaries
%for the GCIP layers.

%search for all manual segmentations in the directory
allSegs = subdir(fullfile(dataDir,'*seg*nii*'));

%process each segmentation
for n = 1:length(allSegs)
    currSeg = allSegs(n).name;
    [currPath,currName,currExt]=fileparts(currSeg);
    
    %check if already processed
    if(isfile(fullfile(currPath,'ManualBoundaries.mat')))
        disp([fullfile(currPath,'ManualBoundaries.mat') ' found, skipping.'])
  %      continue;
    end
    
    %load the segmentation .nii
    segImg = niftiread(currSeg);
    segImg = rot90(segImg,-2);
    segImg = permute(segImg,[2 1 3]);
    
    %load the original montaged .nii
    OCTImg = niftiread(strrep(strrep(currSeg,'_Segmentation',''),'.gz',''));
    OCTImg = rot90(OCTImg,-2);
    OCTImg = permute(OCTImg,[2 1 3]);
    
    %set the three pieces of the montage
    im1_ = OCTImg(:,:,1);
    im2_ = OCTImg(:,:,2);
    im3_ = OCTImg(:,:,3);
    
    %calculate the overlapping regions and get rid of nans
    overlap= ((im1_ ~= 0) & (~isnan(im1_))) + ((im2_ ~= 0) & (~isnan(im2_))) + ((im3_ ~= 0) & (~isnan(im3_)));
    overlap(overlap == 0) = 1;
    
    im1_(isnan(im1_))=0;
    im2_(isnan(im2_))=0;
    im3_(isnan(im3_))=0;
    
    %calculate the average montage for the orginal image
    octImAvg = (im1_+im2_+im3_)./(overlap);
    %cacluate the segmentation montage
    SegAvg = max(segImg,[],3);
    
    %size of the montage
    YN = size(SegAvg,1);
    XN = size(SegAvg,2);
    
    %calculate binary mask and find connected components 
    SegAvg_bw = SegAvg > 0;
    CC=bwconncomp(SegAvg_bw);
    
    %we expect exactly 3 connect components separated by the fovea and
    %optic disk, return error if more or fewer pieces are found
    if(CC.NumObjects ~= 3)
        disp(['Error:' fullfile(currPath,'ManualBoundaries.mat') ' incorrect number of segments, skipping.'])
        continue;
    end
    
    %now we calculate the boundaries
    boundaries = [];    
    overlay = cat(3,octImAvg,octImAvg,octImAvg);%manual segmentation overlay
    overlaySmooth = cat(3,octImAvg,octImAvg,octImAvg);%smoothed segmentation overlay

    %calculate for each piece c
    for c = 1:3
        currSect = zeros(YN,XN);
        currSect(CC.PixelIdxList{c}) = SegAvg(CC.PixelIdxList{c});
        
        %these are the boundaries
        Top=[];%Inner GC boundary
        Mid=[];%GC_IP boundary
        Bot=[];%Outer IP boundary
        
        %for each Ascan we expect to find one of these boundaries
        for x = 1:XN
            ty=find(currSect(:,x)>0,1,'first');
            my=find(currSect(:,x)>1,1,'first');
            by=find(currSect(:,x)>0,1,'last');
            
            %save boundary location (X,Y) at each Ascan
            if(~isempty(ty))
                Top = [Top ; x ty];
            end
            if(~isempty(my))
                Mid = [Mid ; x my];
            end
            
            if(~isempty(by))
                Bot = [Bot ; x by];
            end
        end
        
        %Pad ends
        padSize = 50;
        TopPadded = padEnds(Top,padSize);
        MidPadded = padEnds(Mid,padSize);
        BotPadded = padEnds(Bot,padSize);
        %fit and Calculate a smoothed version of tehse boundaries
        PolyLevel = 15;
        [TopP, TopS, TopMu] = polyfit(TopPadded(:,1),TopPadded(:,2),PolyLevel);
        [MidP, MidS, MidMu] = polyfit(MidPadded(:,1),MidPadded(:,2),PolyLevel);
        [BotP, BotS, BotMu] = polyfit(BotPadded(:,1),BotPadded(:,2),PolyLevel);
        
        TopYSmooth= polyval(TopP,Top(:,1),TopS, TopMu);
        MidYSmooth= polyval(MidP,Mid(:,1),MidS, MidMu);
        BotYSmooth= polyval(BotP,Bot(:,1),BotS, BotMu);
        
        %this is our smoothed data
        TopSmooth = [Top(:,1) TopYSmooth];
        MidSmooth = [Mid(:,1) MidYSmooth];
        BotSmooth = [Bot(:,1) BotYSmooth];
        
        %draw the boundaries as overlay on the OCT montage
        for t= 1:length(Top)
            overlay(Top(t,2),Top(t,1),:) = [0 1 0]';
            %sometimes smoothing extends past the image, so we ignore these
            %case
            if(round(TopSmooth(t,2)) <=0 || round(TopSmooth(t,2)) >YN)
                continue
            end
            overlaySmooth(round(TopSmooth(t,2)),TopSmooth(t,1),:) = [0 1 0]';
        end
        %repeat for the other boundaries
        for t= 1:length(Mid)
            overlay(Mid(t,2),Mid(t,1),:) = [1 1 0]';
            if(round(MidSmooth(t,2)) <=0 || round(MidSmooth(t,2)) >YN)
                continue
            end
            overlaySmooth(round(MidSmooth(t,2)),MidSmooth(t,1),:) = [1 1 0]';
            
        end
        
        for t= 1:length(Bot)
            overlay(Bot(t,2),Bot(t,1),:) = [1 0 0]';
            if(round(BotSmooth(t,2)) <=0 || round(BotSmooth(t,2)) >YN)
                continue
            end
            overlaySmooth(round(BotSmooth(t,2)),BotSmooth(t,1),:) = [1 0 0]';
            
        end
        
        %this is our output data structure saving all three boundaries and
        %for all 3 sections in each image
        boundaries(c).GC_Inner = Top;
        boundaries(c).GC_IP = Mid;
        boundaries(c).IP_Outer = Bot;
        
        boundariesSmooth(c).GC_Inner = TopSmooth;
        boundariesSmooth(c).GC_IP = MidSmooth;
        boundariesSmooth(c).IP_Outer = BotSmooth;
        
    end
    %save the overlay images and the boundary data
    imwrite(overlay,fullfile(currPath,'AvgAll_Trans_wManual.tif'),'tif')
    imwrite(overlaySmooth,fullfile(currPath,['AvgAll_Trans_wManualSmooth_Padded_p' num2str(PolyLevel) '.tif']),'tif')
    save(fullfile(currPath,'ManualBoundaries.mat'),'boundaries','boundariesSmooth','currSeg')
end