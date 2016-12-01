% #############################################################
% Script to reproduce Figure 1 from Turpin/Mckendrick 2015.
% 
% The goal is to produce a plot of the cumulative retinal
% ganglion cell and receptive field counts as a function of
% eccentricity (mm). 
% 
% Everthing should be in units of mm, mm^2, or cell count
% Conversion between units found at the end of Watson 2014
% and fit from Drasdo 1974.
% 
% MAB 2016
%
% Notes:
% Curcio and Allen 1990 reference frames is the retinal feild
%       Data loaded from .xml is in the retinal coordinate frame
%       No extra steps neded to reference the data.
% Drasdo 2007 usues the visual field refrerence frame
%       Indexing the output of densityRF():
%       ex. RFdensity = densityRf(radDeg,smpPerDeg,'OFF');
%            midPoint = center pixel of output image
%            Temporal retina = left of midPoint
%            Nasal Retina    = right of midPoint
%            Supior Retina   = above midPoint
%            Inferior Retina = below midPoint
% Turpin McKendrick 2015 reference frame is the visual field
%
% Watson 2014 reference frame is the visual field
%
% THIS CODE SETS THE REFERENCE FRAME TO THE RETINAL FIELD
%
% #############################################################

%% Generate the Receptive Field Density From the Drasdo 2007:
% Set Parameters 
radDeg      = 20; % Radius in Degrees
smpPerDeg   = 2; % Samples per Degree
radMM = 5;
smpPerMM = 6;
sectorAngle = 6;
rotDegs =0:sectorAngle:360-sectorAngle;

%% Presize output matrix
out_smps = 1/smpPerMM; % parameter from Turpin code
out_radii = 0:out_smps:radMM; % vector of radii to match Turpin code
mapFull = nan(2*length(out_radii)-1,2*length(out_radii)-1,length(rotDegs));
mapFullMidPt = round(length(mapFull)/2);

% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
[RGCdensity,sampleBase_RGC_mm]= densityRGC(radMM,smpPerMM,'OFF');

[RFdensity,sampleBase_RF_deg] = densityRf(radDeg,smpPerDeg,'OFF'); % Generates a 2D Receptive Field Density plot 


for i = 1:length(rotDegs)
    
    [RGCdenisty_mmSq,RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDegs(i));
    
    % Convert the RF counts (and sample base) from degress / deg^2 to mm and cells/mm^2
    RFdensity_mm = convert_degSq_to_mmSq(sampleBase_RF_deg, RFdensity_sqDeg);
    sampleBase_RF_mm=convert_deg_to_mm(sampleBase_RF_deg);
    
    
    mapFull(mapFullMidPt,mapFullMidPt:end,i) = calcDisplacement(RFdensity_mm,sampleBase_RF_mm,RGCdenisty_mmSq,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle,'OFF');
    
    mapFull(:,:,i) = imrotate(mapFull(:,:,i),-1.*rotDegs(i),'crop','nearest');
    
end

dispMap = nanmean(mapFull,3);
%figure;plot(-radMM:1/smpPerMM:radMM,dispMap(31,:),'r');xlabel('eccentricity (mm)'),ylabel('displacement (mm)')