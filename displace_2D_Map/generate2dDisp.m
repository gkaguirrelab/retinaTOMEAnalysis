function dispMap = generate2dDisp(RFdensity,RGCdensity,sampleBase_RF_deg,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle,plot)

rotDegs =0:sectorAngle:360-sectorAngle;

%% Presize output matrix
out_smps = 1/smpPerMM; % parameter from Turpin code
out_radii = 0:out_smps:radMM; % vector of radii to match Turpin code
mapFull = nan(2*length(out_radii)-1,2*length(out_radii)-1,length(rotDegs));
mapFullMidPt = round(length(mapFull)/2);



%% Generate a displacemnt vector for each sector in the retinal feild (nummber of sectors to 360/sectorAngle)
 % This will be stored in a 3D martix with (:,:,i) slice containing a single displacemnt vector 
for i = 1:length(rotDegs)
     % Rotate the RFdensity and the RGCdensity by theta = i*sectorAngle to extract the radial densitiy data    
    [RGCdenisty_mmSq,RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDegs(i));
    
    % Convert the RF counts (and sample base) from degress / deg^2 to mm and cells/mm^2
    RFdensity_mm = convert_degSq_to_mmSq(sampleBase_RF_deg, RFdensity_sqDeg);
    sampleBase_RF_mm=convert_deg_to_mm(sampleBase_RF_deg);
    
    % Calculate the displacment 
    mapFull(mapFullMidPt,mapFullMidPt:end,i) = calcDisplacement(RFdensity_mm,sampleBase_RF_mm,RGCdenisty_mmSq,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle,'OFF');
    
    % radially arrange the displacemnt vectors across the 3rd dim
    mapFull(:,:,i) = imrotate(mapFull(:,:,i),-1.*rotDegs(i),'crop','nearest');
    
end

%% Create 2D map series of radial images
nanDispMap = nanmean(mapFull,3);

%% fill in the nans
dispMap = fillNansInMap(nanDispMap,radMM,smpPerMM);

if strcmp(plot,'full')
    figure;plot(-radMM:1/smpPerMM:radMM,dispMap(31,:),'r');xlabel('eccentricity (mm)'),ylabel('displacement (mm)')

    figure;plot(convert_mm_to_deg(-radMM:1/smpPerMM:radMM),convert_mm_to_deg(dispMap(31,:)),'r');xlabel('eccentricity (deg)'),ylabel('displacement (deg)')
end

end
