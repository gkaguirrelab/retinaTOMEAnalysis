function calcExtOCTGCIPThickness(dataDir)
%Calculate the thickness of the GC and IP layers for the extended OCT 
%with respect to a constant X support
%The function produces .mat file in the same directory as the data with the
%following variables:
% XPos_Pixels[1xP] - The location of the measurements (which is now the same for all subjects, so a single vector)
% XPos_Degs[1xP] - Same as above, but in units of degrees
% GCthicknessValuesAtXPos_um[Nx2xP] - GC thickness values (Dims: 1-Subject, 2-OD/OS, 3-Layer Thickness Values)
% IPthicknessValuesAtXPos_um[Nx2xP] - IP thickness values (Dims: 1-Subject, 2-OD/OS, 3-Layer Thickness Values)
%
% Other relevant variables:
% Sub_AScanResolution_um[Nx2] - resolution (in micron) of the ascan for each subject's image (by eye)
% dataAvailable[Nx2] - flag of whether the data is available for the subject (DIM1) and eye (DIM2, {OD,OS}).
% degPerPixel[1] - static conversion ratio between deg/pixel, according to the machine.
% subIDs[Nx1] - subject IDs, ordered the same as the first dimension of the thickness variables
allFiles = dir(fullfile(dataDir,'1*'));

N = length(allFiles);

degPerPixel=30/1536;%30 degree fov
XPos_Pixels = -1280:1280;%this covers a 25 degrees around the fovea
XPos_Degs = XPos_Pixels*degPerPixel;%same locations in degrees
GCthicknessValuesAtXPos_um = zeros(N,2,length(XPos_Pixels));%(SUB,{OD/OS},THICKNESS)
IPthicknessValuesAtXPos_um = zeros(N,2,length(XPos_Pixels));%(SUB,{OD/OS},THICKNESS)
Sub_AScanResolution_um = zeros(N,2);%Resolution
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
        
        %Save in Microns
        Sub_AScanResolution_um(n,s) = AscanRes*1000;
        GCthicknessValuesAtXPos_um(n,s,:) = minthicknessGC*AscanRes*1000;
        IPthicknessValuesAtXPos_um(n,s,:) = minthicknessIP*AscanRes*1000;

    end
end

save(fullfile(dataDir,"GCIP_thicknessesByDeg.mat"),'subIDs','dataAvailable','degPerPixel','eyeSideIndex','XPos_Pixels','XPos_Degs','GCthicknessValuesAtXPos_um','IPthicknessValuesAtXPos_um','Sub_AScanResolution_um')
