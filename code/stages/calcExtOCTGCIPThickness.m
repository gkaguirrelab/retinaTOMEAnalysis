function calcExtOCTGCIPThickness(dataDir)
%calculate the thickness of the GC and IP layers for the extended OCT 
%with respect to a constant X support
allFiles = dir(fullfile(dataDir,'1*'));

N = length(allFiles);

degPerPixel=30/1536;%30 degree fov
XPos_Pixels = -1280:1280;%this covers a 25 degrees around the fovea
XPos_Degs = XPos_Pixels*degPerPixel;%same locations in degrees
GCthicknessValuesAtXPos = zeros(N,2,length(XPos_Pixels));%(SUB,{OD/OS},THICKNESS)
IPthicknessValuesAtXPos = zeros(N,2,length(XPos_Pixels));%(SUB,{OD/OS},THICKNESS)
eyeSideIndex = {'OD', 'OS'};
dataAvailable = zeros(N,2);%(SUB, {OD/OS})
subIDs = char({allFiles(:).name});
for n = 1:N
    
    for s =1:2
        if (s == 1)
            %OD
            BoundaryFile = fullfile(dataDir,allFiles(n).name,'OD','ManualBoundaries.mat');
            segTemporal = 1;
            segCenter = 2;
            segNasal = 3;
        else
            BoundaryFile = fullfile(dataDir,allFiles(n).name,'OS','ManualBoundaries.mat');
            segTemporal = 2;
            segCenter = 3;
            segNasal = 1;
        end

        if(~isfile(BoundaryFile))
            continue
        else
            dataAvailable(n,s) = 1;
        end
        load(BoundaryFile);

        
        %use segment edges to find fovea/disc centers
        [foveaL_Loc_x, foveaL_Loc_ind] = max(boundariesSmooth(segTemporal).GC_IP(:,1));
        foveaL_Loc_y = boundariesSmooth(segTemporal).GC_IP(foveaL_Loc_ind,2);
        [foveaR_Loc_x, foveaR_Loc_ind] = min(boundariesSmooth(segCenter).GC_IP(:,1));
        foveaR_Loc_y = boundariesSmooth(segCenter).GC_IP(foveaR_Loc_ind,2);
        
        [discL_Loc_x, discL_Loc_ind] = max(boundariesSmooth(segCenter).GC_IP(:,1));
        discL_Loc_y = boundariesSmooth(segCenter).GC_IP(discL_Loc_ind,2);
        [discR_Loc_x, discR_Loc_ind] = min(boundariesSmooth(segNasal).GC_IP(:,1));
        discR_Loc_y = boundariesSmooth(segNasal).GC_IP(discR_Loc_ind,2);
        
        
        %find the fovea location
        foveaLoc_x = mean([foveaL_Loc_x foveaR_Loc_x]);
        foveaLoc_y = mean([foveaL_Loc_y foveaR_Loc_y]);
        
        
        %calculate mask where we actually have boundary information
        support_mask = zeros(size(XPos_Pixels));
        GC_InnerAll=[];
        GC_IPAll = [];
        IP_OuterAll =[];
        
        for c=1:3
            %center boundaries points
            GC_InnerSeg_c = boundariesSmooth(c).GC_Inner - repmat([foveaLoc_x foveaLoc_y],size(boundariesSmooth(c).GC_Inner,1),1);
            GC_IPSeg_c = boundariesSmooth(c).GC_IP - repmat([foveaLoc_x foveaLoc_y],size(boundariesSmooth(c).GC_IP,1),1);
            IP_OuterSeg_c = boundariesSmooth(c).IP_Outer - repmat([foveaLoc_x foveaLoc_y],size(boundariesSmooth(c).IP_Outer,1),1);
            
            segXAll = [GC_InnerSeg_c(:,1); GC_IPSeg_c(:,1); IP_OuterSeg_c(:,1)];
            segXmin = min(segXAll);
            segXmax = max(segXAll);
            
            I = find(XPos_Pixels>=segXmin & XPos_Pixels<=segXmax);
            support_mask(I) = 1;
            
            %Group 3 pieces together
            GC_InnerAll=[GC_InnerAll;GC_InnerSeg_c];
            GC_IPAll = [GC_IPAll;GC_IPSeg_c];
            IP_OuterAll =[IP_OuterAll;IP_OuterSeg_c];
        end
        
        %interplote boundary locations to support location
        GC_Inner_YAtXPos = interp1(GC_InnerAll(:,1),GC_InnerAll(:,2),XPos_Pixels);
        GC_IP_YAtXPos = interp1(GC_IPAll(:,1),GC_IPAll(:,2),XPos_Pixels);
        IP_Outer_YAtXPos = interp1(IP_OuterAll(:,1),IP_OuterAll(:,2),XPos_Pixels);
        
        %remove locations with no data
        GC_Inner_YAtXPos(~support_mask) = 0;
        GC_IP_YAtXPos(~support_mask) = 0;
        IP_Outer_YAtXPos(~support_mask) = 0;
        
        GC_InnerFinal = [XPos_Pixels', GC_Inner_YAtXPos'];
        GC_IPFinal = [XPos_Pixels', GC_IP_YAtXPos'];
        IP_OuterFinal = [XPos_Pixels', IP_Outer_YAtXPos'];
        
        %calculate thickness
        minthicknessGC = boundaryDist(GC_IPFinal,GC_InnerFinal,1);
        minthicknessIP = boundaryDist(GC_IPFinal,IP_OuterFinal,1);
        
        GCthicknessValuesAtXPos(n,s,:) = minthicknessGC;
        IPthicknessValuesAtXPos(n,s,:) = minthicknessIP;

    end
end

save(fullfile(dataDir,"GCIP_thicknessesByDeg_7_18_2019.mat"),'subIDs','dataAvailable','degPerPixel','eyeSideIndex','XPos_Pixels','XPos_Degs','GCthicknessValuesAtXPos','IPthicknessValuesAtXPos')
