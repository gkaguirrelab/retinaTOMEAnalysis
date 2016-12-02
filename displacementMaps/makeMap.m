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


[RGCdensity,sampleBase_RGC_mm]= densityRGC(radMM,smpPerMM,'OFF');

% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
[RFdensity,sampleBase_RF_deg] = densityRf(radDeg,smpPerDeg,'OFF'); % Generates a 2D Receptive Field Density plot 

dispMap = generate2dDisp(RFdensity,RGCdensity,sampleBase_RF_deg,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle,'full');
end