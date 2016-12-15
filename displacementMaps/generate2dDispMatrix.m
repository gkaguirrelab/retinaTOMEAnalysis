function dispMap = generate2dDispMatrix(RFdensity,RGCdensity,sampleBase_RF_deg,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle)

rotDegs =0:sectorAngle:360-sectorAngle;

%% Presize output matrix
out_smps = 1/smpPerMM; % parameter from Turpin code
out_radii = 0:out_smps:radMM; % vector of radii to match Turpin code
dispMap = nan(length(rotDegs),length(out_radii));




%% Generate a displacemnt vector for each sector in the retinal feild (nummber of sectors to 360/sectorAngle)
 % This will be stored in a 3D martix with (:,:,i) slice containing a single displacemnt vector 
for i = 1:length(rotDegs)
     % Rotate the RFdensity and the RGCdensity by theta = i*sectorAngle to extract the radial densitiy data    
    [RGCdensity_mmSq,RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDegs(i));
    
    % Convert the RF counts (and sample base) from degress / deg^2 to mm and cells/mm^2
    RFdensity_mm = convert_degSq_to_mmSq(sampleBase_RF_deg, RFdensity_sqDeg);
    sampleBase_RF_mm=convert_deg_to_mm(sampleBase_RF_deg);
    
    % Calculate the displacment 
    dispMap(i,:) = calcDisplacement(RFdensity_mm,sampleBase_RF_mm,RGCdensity_mmSq,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle,'OFF');
   
end

end
