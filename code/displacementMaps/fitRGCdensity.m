function [RGCdensity,sampleBase_RGC_mm]= fitRGCdensity(radMM,smpPerMM,interp,verbose)
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
polyFunctionNasal = fit(log(ecc_mm(2:end)),log(rgcDenisty_mmSq_nasal(2:end)),'poly8','Exclude', find(isnan(rgcDenisty_mmSq_nasal(2:end))));
p = polyfit(log(ecc_mm(2:end)),log(rgcDenisty_mmSq_nasal(2:end)),14)


polyFunctionSuperior = fit(ecc_mm,rgcDensity_mmSq_superior,'smoothingspline', 'Exclude',find(isnan(rgcDensity_mmSq_superior)),'SmoothingParam', 1);
polyFunctionTemporal = fit(ecc_mm,rgcDensity_mmSq_temporal,'smoothingspline', 'Exclude',find(isnan(rgcDensity_mmSq_temporal)),'SmoothingParam', 1);
polyFunctionInferior = fit(ecc_mm,rgcDenisty_mmSq_inferior,'smoothingspline', 'Exclude',find(isnan(rgcDenisty_mmSq_inferior)),'SmoothingParam', 1);

figure;
subplot(2,1,1)
hold on
plot(log(ecc_mm(2:end)),log(rgcDenisty_mmSq_nasal(2:end)));
plot(log(ecc_mm(2:end)),polyFunctionNasal(log(ecc_mm(2:end))))
errorbar(x,y,err,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')


subplot(2,1,2)
figure
hold on
%plot(ecc_mm(2:end),rgcDenisty_mmSq_nasal(2:end));
plot(ecc_mm(2:end),exp(polyFunctionNasal(log(ecc_mm(2:end)))),'r')
errorbar(ecc_mm(2:end),rgcDenisty_mmSq_nasal(2:end),rgcDenisty_SD_nasal(2:end),'-s','MarkerSize',3,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue')







