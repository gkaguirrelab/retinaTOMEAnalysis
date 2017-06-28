function [dispMap]  = makeMap(varargin)
% function [dispMap]  = makeMap(varargin)
%
% Description:
%  Generate a two-dimensional map of the displacement of retinal ganglion
%  cells from their receptive field centers.
%
%  The solution is provided as millimeters of displacement towards the
%  fovea as a function of eccentrict distance from the fovea in mm. This is
%  calculated for each of many radial directions.
%
% Inputs (given as key/value pairs):
%
%   radMM - The radius (in mm) that defines the spatial support of the
%     displacement map
%   smpPerMM - The number of samples per mm at which to perform the
%     displacement calculation along the radial.
%   sectorAngle - The spacing in degrees of each radial.
%
% Outputs:
%
%   dispMap - A m x n matrix, where m is the number of radials (360 /
%     sectorAngle) and n is the number of samples along each radial (radMM
%     * smpPerMM +1). The value given at each point is the displacement in
%     mm of the RGC cell bodies back towards the fovea along that radial
%     angle to align the cell bodies with their receptive fields.
%
% Notes:
%
%   The radials are ordered in increasing theta values, starting from zero,
%   and reaching a maximum of [360 - sectorAngle]. Theta of zero
%   corresponds to horizontal radial along the nasal retina. A theta of 90
%   corresponds to the vertical radial along the superior retina.
%
% MAB 2016

%% Parse vargin for options passed here
p = inputParser;

p.addParameter('radMM',5,@isnumeric);
p.addParameter('smpPerMM',2,@isnumeric);
p.addParameter('sectorAngle',3,@isnumeric);

p.parse(varargin{:});

% Create a map of RGC density within the retinal extent defined by the
% input parameters.
[RGCdensity,sampleBase_RGC_mm]= densityRGC(p.Results.radMM,p.Results.smpPerMM,'spline','OFF');

% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
[RFdensity,sampleBase_RF_deg] = densityRf(convert_mm_to_deg(p.Results.radMM),((p.Results.radMM*p.Results.smpPerMM)/convert_mm_to_deg(p.Results.radMM)),'spline','OFF'); % Generates a 2D Receptive Field Density plot 

dispMap = generate2dDispMatrix(RFdensity,RGCdensity,sampleBase_RF_deg,sampleBase_RGC_mm,p.Results.radMM,p.Results.smpPerMM,p.Results.sectorAngle);
end