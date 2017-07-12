function [ecc_mm,outParams,RGCdensityFit, scaleData] = fitRGCdensityDev(angle)
% fitRGCdensity -- Estimates a RGC cell body density function at a given angle on the retina.
%
% Description:
%   This function returns a function that esimates the retinal ganglion cell
%   cell body density as a function of eccentricity in mm .This uses the 
%   Curcio and Allen (1990) data that measured cell density along the 4 meridians.
%   This returns a function that estimates this density at the desired
%   input angle by taking a weighted average of the fit parameters.
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

%% Load the RGC Density Data from Curcio and Allen 1990:
% Curcio and Allen obtained measurements of the denisty of all RGC classes
% within 6 human retinas at a set of positions relative to the fovea. These
% data were provided online in 2013.
curcio_data = fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio/curcio_4meridian.mat']);
data=load(curcio_data);
data=data.data;

% Here we distribute the contents of the Curcio measurements across some
% variables:
%  ecc_mm - defines the distance in mm from the fovea along each of the
%           radials for which we have denisty measurements.
%  temp / sup / nasal / inferior - RGC density (counts / mm^2) at each
%           position along each radial. The labels correspond to retinal
%           coordinate frame (i.e., nasal refers to the nasal retina).
%           This can be verified by observing that there is an interruption
%           in the count data for the nasal meridian corresponding to the
%           blind spot.
ecc_mm = data(:,1);

rgcDensity_mmSq_temporal = data(:,2);
rgcDensity_mmSq_superior = data(:,4);
rgcDensity_mmSq_nasal = data(:,6);
rgcDensity_mmSq_inferior = data(:,8);

rgcDensity_SD_temporal = data(:,3);
rgcDensity_SD_superior = data(:,5);
rgcDensity_SD_nasal = data(:,7);
rgcDensity_SD_inferior = data(:,9);

% Adjust the RGC counts to reflect the fraction that is midget cells
midgetFraction = midgetFractionByEccen(convert_mm_to_deg(ecc_mm));
midget_rgcDensity_mmSq_superior = rgcDensity_mmSq_superior.* midgetFraction;

% set the scale to normalize the data so that the  Frechet PDF fits well.
% ## THIS MUST ALSO BE APPLIED TO THE RF DENSITY DATA ##
%scaleData = 2*max([rgcDensity_mmSq_temporal;midget_rgcDensity_mmSq_superior;rgcDensity_mmSq_nasal;rgcDensity_mmSq_inferior]);
scaleData = 2*max([midget_rgcDensity_mmSq_superior]);
% Normalize the data
norm_rgcDensity_temporal = rgcDensity_mmSq_temporal./scaleData;
norm_rgcDensity_superior = midget_rgcDensity_mmSq_superior./scaleData;
norm_rgcDensity_nasal = rgcDensity_mmSq_nasal./scaleData;
norm_rgcDensity_inferior = rgcDensity_mmSq_inferior./scaleData;


% Fit a Frechet to each of the cardinal radial directions across
% the sampled RGC densities. The result is a set of function handles that
% relate continuous distance (in mm) from the fovea to RGC density (in
% counts / mm^2).

%% THERE IS A BUG IN HOW THE ANGLE IS BEING RELATED TO THE MIXING OF THE MERIDIA

%% REPLACE WITH SWITCH/CASE cause it's prettier

% Take a weighed average of the parameters of the fit polynomial. the
% weights are the fraction of the input angle for both meridians that flank the
% angle of interest. 
if angle >= 0 && angle < 90
    nasalFrac = angle/90;
    superiorFrac = 1 - nasalFrac;
    weightedAvgData = nasalFrac.* norm_rgcDensity_nasal + superiorFrac.*norm_rgcDensity_superior;
elseif angle >= 90 && angle < 180
    superiorFrac = (angle-90)/90;
    temporalFrac = 1 - superiorFrac;
    weightedAvgData = superiorFrac.*norm_rgcDensity_superior + temporalFrac.*norm_rgcDensity_temporal;
elseif angle >= 180 && angle < 270
    temporalFrac = (angle-180)/90;
    inferiorFrac = 1 - temporalFrac;
    weightedAvgData = temporalFrac.*norm_rgcDensity_temporal + inferiorFrac.*norm_rgcDensity_inferior;
elseif angle >= 270 && angle < 360
    inferiorFrac = (angle-270)/90;
    nasalFrac = 1 - inferiorFrac;
    weightedAvgData = inferiorFrac.*norm_rgcDensity_inferior + nasalFrac.*norm_rgcDensity_nasal;
end

%% SHORT-CIRCUITED TO ONLY RETURN THE SUPERIOR MERIDIAN FIT
%[outParams, RGCdensityFit] = fitFrechetToRGCDensity(ecc_mm, weightedAvgData, ones(size(ecc_mm)));
[outParams, RGCdensityFit] = fitFrechetToRGCDensity(ecc_mm, norm_rgcDensity_superior, ones(size(ecc_mm)));

end


function midgetFraction = midgetFractionByEccen(ecc_mm)
% function midgetFraction = midgetFractionByEccen(ecc_mm)
%
% This function returns the fraction of midget retinal ganglion cells in
% relation to the entire count of RGCs at a given distance (in mm) from the
% fovea. The equation is taken from Watson JoV 2016, 

f0 = 0.8928;
rm = 41.03;
midgetFraction = f0.*(1+(ecc_mm/rm)).^-1;
end % function

