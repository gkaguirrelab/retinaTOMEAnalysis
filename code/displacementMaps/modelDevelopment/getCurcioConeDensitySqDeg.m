function [supportPosDeg,coneDensitySqDeg] = getCurcioConeDensitySqDeg(angle)
% fitConeDensity -- Loads and fits the Curcio 1990 cone density data
%
% Description:
%   This function returns a function that fits the cone density data 
%   reported in:
%
%      Curcio et al (1990) Human photoreceptor topography. J Comp Neurology
%
%
% Inputs:
%   angle      = The dedsired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
%   polyOrder  = The order of the polynomial used for fitting.
% Outputs:
%   expFitOut = Function that estimates the RGC cell body denstity at 
%               the input anlge as a funciton of eccentricity (mm).
%
% MAB 2017

%% Load the Cone Density Data from Curcio et al 1990:
% Curcio and Allen obtained measurements of the denisty of cones
% within 6 human retinas at a set of positions relative to the fovea.
% Loading this matlab brings the variable "curcioConeDensityPerSqMm" into
% memory
curcioConeDataFile = fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio_1990_JCompNeurol_HumanPhotoreceptorTopography/curcioConeDensityPerSqMm.mat']);
load(curcioConeDataFile);

meridianNames = {'superior','inferior','temporal','nasal'};

% Convert mm to deg, and mm^2 to deg^2
curcioConeDensityPerSqDeg.support = convert_mm_to_deg(curcioConeDensityPerSqMm.support);

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
    otherwise
        error('Only have the cardinal meridians right now');
end

end

