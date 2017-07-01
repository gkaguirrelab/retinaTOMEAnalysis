function [ecc_mm,polyFitOut] = fitRGCdensity(angle,polyOrder)
% densityRGC -- Retinal Ganglion Cell 2D density map.
%
% Description:
%   This function produces a two dimensional map of retinal ganglion cell
%   density.This uses the Curcio and Allen (1990) retinal ganglion cell
%   denstiy data which has RGC density measurses (per mm^2) along the four
%   meridians and radially interpolates betwwen the arms of the meridians.
%   meridians data
%
% Inputs:
%   radMM    = desired radius of the map in mm.
%   smpPerMM = how many samples per mm.
%   interp   = interpoation method for interp1 ('spline', 'linear',... see
%              help interp1 for more intput options).
%   verbose  = option to plot density map (1 = plot, 0 = no plot).
%
% Outputs:
%   RGCdensity        = Map of retinal ganglion cell density.
%   sampleBase_RGC_mm = Sample point positions from 0 to radMM in mm.
%
% MAB 2016

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
rgcDenisty_mmSq_nasal = data(:,6);
rgcDenisty_mmSq_inferior = data(:,8);

rgcDensity_SD_temporal = data(:,3);
rgcDensity_SD_superior = data(:,5);
rgcDenisty_SD_nasal = data(:,7);
rgcDenisty_SD_inferior = data(:,9);

% Fit a smoothing spline to each of the cardinal radial directions across
% the sampled RGC densities. The result is a set of function handles that
% relate continuous distance (in mm) from the fovea to RGC density (in
% counts / mm^2).

polyFunctionNasal = fit(log(ecc_mm),log(rgcDenisty_mmSq_nasal),['poly' num2str(polyOrder)],'Exclude', [find(isnan(rgcDenisty_mmSq_nasal)) find(isinf(log(rgcDenisty_mmSq_nasal)))']);
polyFunctionSuperior = fit(log(ecc_mm),log(rgcDensity_mmSq_superior),['poly' num2str(polyOrder)],'Exclude', [find(isnan(rgcDensity_mmSq_superior)) find(isinf(log(rgcDensity_mmSq_superior)))']);
polyFunctionTemporal = fit(log(ecc_mm),log(rgcDensity_mmSq_temporal),['poly' num2str(polyOrder)],'Exclude', [find(isnan(rgcDensity_mmSq_temporal))' find(isinf(log(rgcDensity_mmSq_temporal)))']);
polyFunctionInferior = fit(log(ecc_mm),log(rgcDenisty_mmSq_inferior),['poly' num2str(polyOrder)],'Exclude', [find(isnan(rgcDenisty_mmSq_inferior))' find(isinf(log(rgcDenisty_mmSq_inferior)))']);

polyFitOut = polyFunctionNasal;

if angle >= 0 && angle < 90;
    nasalFrac = angle/90;
    superiorFrac = 1 - nasalFrac;
    for i = 1:polyOrder
        eval(sprintf('polyFitOut.p%s = nasalFrac.*polyFunctionNasal.p%s + superiorFrac.*polyFunctionSuperior.p%s;',num2str(i),num2str(i),num2str(i)))
    end
elseif angle >= 90 && angle < 180;
    superiorFrac = (angle-90)/90;
    temporalFrac = 1 - superiorFrac;
    for i = 1:polyOrder
        eval(sprintf('polyFitOut.p%s = superiorFrac.*polyFunctionSuperior.p%s + temporalFrac.*polyFunctionTemporal.p%s;',num2str(i),num2str(i),num2str(i)))
    end
elseif angle >= 180 && angle < 270;
    temporalFrac = (angle-180)/90;
    inferiorFrac = 1 - temporalFrac;
    for i = 1:polyOrder
        eval(sprintf('polyFitOut.p%s = temporalFrac.*polyFunctionTemporal.p%s + inferiorFrac.*polyFunctionInferior.p%s;',num2str(i),num2str(i),num2str(i)))
    end
elseif angle >= 270 && angle < 360;
    inferiorFrac = (angle-270)/90;
    nasalFrac = 1 - inferiorFrac;
    for i = 1:polyOrder
        eval(sprintf('polyFitOut.p%s = inferiorFrac.*polyFunctionInferior.p%s + nasalFrac.*polyFunctionNasal.p%s;',num2str(i),num2str(i),num2str(i)))
    end
end



end






