function [ rgcDensitySqDeg, supportPosDeg ] = getCurcioRGCDensityByEccen(polarAngle)
% getCurcioMidgetRGCDensityByEccen(angle)
%
% Curcio and Allen obtained measurements of the density of all RGC classes
% within 6 human retinas at a set of positions relative to the fovea. These
% data were provided online in 2013.
%
% Here, we load these data and convert from mm to degrees.
% This routine does not adjust for the midget fraction at each eccentricity
%
% Curcio and Allen 1990 reference frame is the retinal feild
%       Data loaded from .xml is in the retinal coordinate frame
%       No extra steps neded to reference the data.
%
% Inputs
%   polarAngle - the meridian to evaluate. Acceptable values are: 
%       (0=nasal; 90=superior; 180=temporal; 270=inferior)
%
% Outputs
%   rgcDensitySqDeg - the density (counts per square
%       degree) of RGCs at each of the positions
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which the RGC density is defined

% Check the input
if sum([0 90 180 270]==polarAngle) ~= 1
    error('The Curcio and Allen data are defined only for the cardinal meridia');
end

% Load the RGC Density Data from Curcio and Allen 1990:
curcioDataFileName = ...
    fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio_1990_JCompNeurol_GanglionCellTopography/curcio_4meridian.mat']);
dataLoad=load(curcioDataFileName);
curcioRGCdensity_mm = dataLoad.data;

% The first column of the data  defines the distance in mm from the fovea
%  along each of the radials for which we have density measurements. Obtain
%  this and convert to degrees
supportPosMm  = curcioRGCdensity_mm(:,1);
supportPosDeg = convert_mm_to_deg(supportPosMm);

% The other columns contain the mean (across eyes) measurements. Use a
% switch statement to grab the column corresponding to the passed angle
switch polarAngle
    case 0 % nasal
        rgcDensitySqDeg = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,6));
        
    case 90 % superior
        rgcDensitySqDeg = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,4));
        
    case 180 % temporal
        rgcDensitySqDeg = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,2));
        
    case 270 % inferior
        rgcDensitySqDeg = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,8));
        
end


end

