function [ midget_rgcDensity_degSq, supportPosDeg ] = getCurcioMidgetRGCDensityByEccen(angle)
% getCurcioMidgetRGCDensityByEccen(angle)
%
% Curcio and Allen obtained measurements of the denisty of all RGC classes
% within 6 human retinas at a set of positions relative to the fovea. These
% data were provided online in 2013.
%
% Here, we load these data, convert from mm to degrees, and then adjust for
% the midget fraction at each eccentricity
%
% Curcio and Allen 1990 reference frame is the retinal feild
%       Data loaded from .xml is in the retinal coordinate frame
%       No extra steps neded to reference the data.
%
% Inputs
%   angle - the meridian to evaluate. Acceptable values are: 
%       (0=nasal; 90=superior; 180=temporal; 270=inferior)
%
% Outputs
%   midget_rgcDensity_degSq - the density (counts per square
%       degree) of midget RGCs at each of the positions
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which the midget RGC density is defined

% Check the input
if sum([0 90 180 270]==angle) ~= 1
    error('The Curcio and Allen data are defined only for the cardinal meridia');
end

% Load the RGC Density Data from Curcio and Allen 1990:
curcioDataFileName = fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio/curcio_4meridian.mat']);
dataLoad=load(curcioDataFileName);
curcioRGCdensity_mm = dataLoad.data;

% The first column of the data  defines the distance in mm from the fovea
%  along each of the radials for which we have density measurements. Obtain
%  this and convert to degrees
supportPosMm  = curcioRGCdensity_mm(:,1);
supportPosDeg = convert_mm_to_deg(supportPosMm);

% The other columns contain the mean (across eyes) measurements. Use a
% switch statement to grab the column corresponding to the passed angle
switch angle
    case 0 % nasal
        rgcDensity_degSq = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,6));
        
    case 90 % superior
        rgcDensity_degSq = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,4));
        
    case 180 % temporal
        rgcDensity_degSq = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,2));
        
    case 270 % inferior
        rgcDensity_degSq = convert_mmSq_to_degSq(supportPosDeg,curcioRGCdensity_mm(:,8));
        
end

% Adjust the RGC counts to reflect the fraction that is midget cells
midgetFraction = midgetFractionByEccen(supportPosDeg);
midget_rgcDensity_degSq = rgcDensity_degSq .* midgetFraction;


end
