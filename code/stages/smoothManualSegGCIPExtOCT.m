function smoothManualSegGCIPExtOCT(dataDir,origFileNameSuffix)
%Purpose: Find all manual segmentations in dataDir with Horizontal Extended
%OCT data and calculate and smooth the the manual segmentation boundaries
%for the GCIP layers.


% Sometimes stray pixels appear in the manual segmentation. If any
% connected cluster of pixels is less than this threshold, remove these
% specks from the segmentation.
componentPixelThresh = 10;

%search for all manual segmentations in the directory
allSegs = subdir(fullfile(dataDir,'*seg*nii'));

%process each segmentation
for n = 1:length(allSegs)
    ValidFlag = 1;
    currSeg = allSegs(n).name;
    [currPath,currName,currExt]=fileparts(currSeg);
    
    %load the segmentation .nii
    
    segImg = niftiread(currSeg);
    segImg = rot90(segImg,-2);
    segImg = permute(segImg,[2 1 3]);
    
    %load the original montaged .nii
    OCTFile = strrep(strrep(currSeg,'_Segmentation',origFileNameSuffix),'.gz','');
    if contains(OCTFile,'11094_OD')
        foo =1 ;
    end
    OCTInfo=niftiinfo(OCTFile);
    OCTImg = niftiread(OCTFile);
    if ndims(OCTImg) == 3
        OCTImg = rot90(OCTImg,-2);
        OCTImg = permute(OCTImg,[2 1 3]);
    else
        OCTImg = flipud(rot90(OCTImg,-1));
    end
    
    %set the three pieces of the montage
    if ndims(OCTImg) == 3
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
        
    else
        
        octImAvg = OCTImg;
        SegAvg = segImg;
        
    end
    
    
    %size of the montage
    YN = size(SegAvg,1);
    XN = size(SegAvg,2);
    
    %calculate binary mask and find connected components
    SegAvg_bw = SegAvg > 0;
    
    % Remove any stray components that are below the size threshold
    CC=bwconncomp(SegAvg_bw);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    if any(numPixels < componentPixelThresh)
        removeIdx = find(numPixels < componentPixelThresh);
        for bb = 1:length(removeIdx)
            SegAvg_bw(CC.PixelIdxList{removeIdx(bb)})=0;
        end
        CC=bwconncomp(SegAvg_bw);
    end
    
    % For the vertical line scans, we expect 2 connected components
    % separated by the fovea, and for the horizontal line scan, we expect
    % three segments, as the nasal retina is also broken up by the optic
    % disc.
    
    if CC.NumObjects~=2 && CC.NumObjects~=3
        disp(['Error:' fullfile(currPath,currSeg) ' incorrect number of segments, skipping.'])
        ValidFlag = 0;
        continue;
    end
    
    boundaries = [];
    rgbDims = 3;
    overlay = double(cat(rgbDims,octImAvg,octImAvg,octImAvg));%manual segmentation overlay
    overlaySmooth = double(cat(rgbDims,octImAvg,octImAvg,octImAvg));%smoothed segmentation overlay
    
    boundaries = [];
    boundariesSmooth = [];
    
    %calculate for each piece c
    for c = 1:CC.NumObjects
        currSect = zeros(YN,XN);
        currSect(CC.PixelIdxList{c}) = SegAvg(CC.PixelIdxList{c});
        
        %these are the boundaries
        Top=[];%Inner GC boundary
        Mid=[];%GC_IP boundary
        Bot=[];%Outer IP boundary
        
        %for each Ascan we expect to find one of these boundaries
        for x = 1:XN
            
            if(~isempty(find(currSect(:,x)>2,1)))
                continue;
            end
            
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
        
        
        if(isempty(Top) || isempty(Mid) || isempty(Bot))
            disp(['Error:' fullfile(currPath,currSeg) ' missing boundaries, skipping.'])
            ValidFlag=0;
            continue;
        end
        
        
        TopPadded = padEnds(Top,padSize);
        MidPadded = padEnds(Mid,padSize);
        BotPadded = padEnds(Bot,padSize);
        %fit and Calculate a smoothed version of these boundaries
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
%         for t= 1:length(Top)
%             overlay(Top(t,2),Top(t,1),:) = [0 1 0]';
%             %sometimes smoothing extends past the image, so we ignore these
%             %case
%             if(round(TopSmooth(t,2)) <=0 || round(TopSmooth(t,2)) >YN)
%                 continue
%             end
%             overlaySmooth(round(TopSmooth(t,2)),TopSmooth(t,1),:) = [0 1 0]';
%             overlaySmooth(round(TopSmooth(t,2))+1,TopSmooth(t,1),:) = [0 1 0]';
%         end
        %repeat for the other boundaries
        for t= 1:length(Mid)
            overlay(Mid(t,2),Mid(t,1),:) = [1 1 0]';
            if(round(MidSmooth(t,2)) <=0 || round(MidSmooth(t,2)) >YN)
                continue
            end
            overlaySmooth(round(MidSmooth(t,2)),MidSmooth(t,1),:) = [1 0 0]';
            overlaySmooth(round(MidSmooth(t,2))+1,MidSmooth(t,1),:) = [1 0 0]';
            
        end
        
        for t= 1:length(Bot)
            overlay(Bot(t,2),Bot(t,1),:) = [1 0 0]';
            
            if(round(BotSmooth(t,2)) <=0 || round(BotSmooth(t,2)) >YN)
                continue
            end
            overlaySmooth(round(BotSmooth(t,2)),BotSmooth(t,1),:) = [1 0 0]';
            overlaySmooth(round(BotSmooth(t,2))+1,BotSmooth(t,1),:) = [1 0 0]';            
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
    AscanRes = OCTInfo.PixelDimensions(2);%resolution of A-scans
    if(ValidFlag) %only save if no errors were found
        save(fullfile(currPath,'ManualBoundaries.mat'),'boundaries','boundariesSmooth','currSeg','AscanRes')
    end
end