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
rotDeg = 15;
% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
[RGCdensity,sampleBase_RGC_mm]= densityRGC(radMM,smpPerMM,'OFF');

RFdensity = densityRf(radDeg,smpPerDeg,'OFF'); % Generates a 2D Receptive Field Density plot 

% Get data from Superior Merdian as a first pass check to validate the
% pipeline of Turpin/McKendrick 
[RGCdenisty_mmSq,RFdensity_sqDeg]= rotAndExrtractMeridian(RFdensity,RGCdensity,rotDeg);

% Calculate a sample base
sampleBase_RF_deg = (0:1/smpPerDeg:radDeg)'; % the eccentricy of each sample of sup_RFdensity in degrees 

% Convert the RF counts (and sample base) from degress / deg^2 to mm and cells/mm^2
RFdensity_mm = convert_degSq_to_mmSq(sampleBase_RF_deg, RFdensity_sqDeg);
sampleBase_RF_mm=convert_deg_to_mm(sampleBase_RF_deg);

Displacement = calcDisplacement(RFdensity_mm,sampleBase_RF_mm,RGCdenisty_mmSq,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle);
