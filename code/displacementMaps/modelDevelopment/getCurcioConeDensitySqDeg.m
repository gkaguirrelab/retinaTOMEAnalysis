function [coneDensitySqDeg, supportPosDeg] = getCurcioConeDensitySqDeg(angle)
% getCurcioConeDensitySqDeg(angle)
%
% This routine returns the cone density data reported in:
%
%      Curcio et al (1990) Human photoreceptor topography. J Comp Neurology
%
% Here, we load these data and convert from mm to degrees.
%
% Inputs:
%   angle      = The dedsired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
% Outputs:
%   coneDensitySqDeg - the density (counts per square
%       degree) of cones at each of the positions
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which the cone density is defined
%

% Check the input
if sum([0 90 180 270]==angle) ~= 1
    error('The Curcio cone data are defined only for the cardinal meridia');
end

%% Load the Cone Density Data from Curcio et al 1990:
% Curcio and Allen obtained measurements of the density of cones
% within 6 human retinas at a set of positions relative to the fovea.
% Loading this matlab brings the variable "curcioConeDensityPerSqMm" into
% memory
curcioConeDataFile = ...
    fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioConeDensityPerSqMm.mat']);
load(curcioConeDataFile);

meridianNames = {'superior','inferior','temporal','nasal'};

% Convert mm to deg, and mm^2 to deg^2
curcioConeDensityPerSqDeg.support = ...
    convert_mm_to_deg(curcioConeDensityPerSqMm.support);

for mm = 1:length(meridianNames)
curcioConeDensityPerSqDeg.(meridianNames{mm}) = ...
    convert_mmSq_to_degSq(curcioConeDensityPerSqDeg.support, curcioConeDensityPerSqMm.(meridianNames{mm}) );
end

supportPosDeg = curcioConeDensityPerSqDeg.support;
switch angle
    case 0
        coneDensitySqDeg = curcioConeDensityPerSqDeg.nasal;
    case 90
        coneDensitySqDeg = curcioConeDensityPerSqDeg.superior;
    case 180
        coneDensitySqDeg = curcioConeDensityPerSqDeg.temporal;
    case 270
        coneDensitySqDeg = curcioConeDensityPerSqDeg.inferior;
end

end % function

