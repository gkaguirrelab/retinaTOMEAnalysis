function dispMap = makeMap(radDeg,smpPerDeg,radMM,smpPerMM,sectorAngle)

% #############################################################
% Script produces a 2D displacement image.
% 
% The goal is to produce an image of the displacement of retinal
% ganglion cell from receptive fields as a function of
% eccentricity (mm) and polar angle. 
% 
% The output is in mm of displacement.
% 
% EXAMPLE INPUT PARAMETERS:
%  radDeg      = 20; % Radius in Degrees
%  smpPerDeg   = 2; % Samples per Degree
%  radMM = 5;
%  smpPerMM = 6;
%  sectorAngle = 6;
%
% Notes:
% 
% THIS CODE SETS THE REFERENCE FRAME TO THE RETINAL FIELD
%
% MAB 2016
% #############################################################

rotDegs =0:sectorAngle:360-sectorAngle;

%% Presize output matrix
out_smps = 1/smpPerMM; % parameter from Turpin code
out_radii = 0:out_smps:radMM; % vector of radii to match Turpin code
mapFull = nan(2*length(out_radii)-1,2*length(out_radii)-1,length(rotDegs));
mapFullMidPt = round(length(mapFull)/2);

% Obtain the ganglion cell density per square mm within the sampling
% area specified (also in mm)
[RGCdensity,sampleBase_RGC_mm]= densityRGC(radMM,smpPerMM,'OFF');

% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
[RFdensity,sampleBase_RF_deg] = densityRf(radDeg,smpPerDeg,'OFF'); % Generates a 2D Receptive Field Density plot 

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

%figure;plot(-radMM:1/smpPerMM:radMM,dispMap(31,:),'r');xlabel('eccentricity (mm)'),ylabel('displacement (mm)')
end