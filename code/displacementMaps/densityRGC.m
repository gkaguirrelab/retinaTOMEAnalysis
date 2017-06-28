function [RGCdensity,sampleBase_RGC_mm]= densityRGC(radMM,smpPerMM,interp,verbose)
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

% Fit a smoothing spline to each of the cardinal radial directions across
% the sampled RGC densities. The result is a set of function handles that
% relate continuous distance (in mm) from the fovea to RGC density (in
% counts / mm^2).
splineFunctionNasal = fit(ecc_mm,rgcDenisty_mmSq_nasal,'smoothingspline','Exclude', find(isnan(rgcDenisty_mmSq_nasal)),'SmoothingParam', 1);
splineFunctionSuperior = fit(ecc_mm,rgcDensity_mmSq_superior,'smoothingspline', 'Exclude',find(isnan(rgcDensity_mmSq_superior)),'SmoothingParam', 1);
splineFunctionTemporal = fit(ecc_mm,rgcDensity_mmSq_temporal,'smoothingspline', 'Exclude',find(isnan(rgcDensity_mmSq_temporal)),'SmoothingParam', 1);
splineFunctionInferior = fit(ecc_mm,rgcDenisty_mmSq_inferior,'smoothingspline', 'Exclude',find(isnan(rgcDenisty_mmSq_inferior)),'SmoothingParam', 1);

% Define the spatial support over which we will interpolate a 2D map of RGC
% density:
%   meridian dim 1 - presize this for output
%   meridian dim 2 - polar angle (in deg)
%   meridian dim 3 - eccentricity (in mm)

meridian = zeros(2*(round(radMM)*round(smpPerMM))+1,2*(round(radMM)*round(smpPerMM))+1);
[~,meridian(:,:,2),meridian(:,:,3)] = createGrid(radMM,smpPerMM);

% interpolate RGC denisty across the spatial support
for i = 1:size(meridian,1)
    for ii = 1:size(meridian,2)
        
        if meridian(i,ii,2) >= 0 && meridian(i,ii,2) <= 90
            VMd = splineFunctionSuperior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = splineFunctionNasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([0,90],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) > 90 && meridian(i,ii,2) <= 180
            VMd = splineFunctionSuperior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = splineFunctionTemporal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([90,180],[VMd,HMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 180 && meridian(i,ii,2) <= 270
            VMd = splineFunctionInferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = splineFunctionTemporal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([180,270],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 270 && meridian(i,ii,2) <=360
            VMd = splineFunctionInferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = splineFunctionNasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([270,360],[VMd,HMd],theta,interp);
            
        end
    end
end 

% NaN points that are beyond the requested radius but are still on the
% Cartesian grid (i.e., the corners of the matrix)
mask=(meridian(:,:,3)<=radMM);
mask=double(mask);
mask(mask == 0) = NaN;
RGCdensity = meridian(:,:,1).*mask;
sampleBase_RGC_mm = 0:1/smpPerMM:radMM;

%% Validate the Output

if strcmp(verbose,'full')
    % 0-Nasal 90-Surperior 180-Temporal 270-Inferior
    rgc = meridian(:,:,1);
    mdPt = round(size(rgc,1)/2);
    xSmpDegPos = 0:1/smpPerMM:radMM;
    xSmpDegNeg = radMM:-1/smpPerMM:0;
    figure;hold on;
    plot(xSmpDegNeg,rgc(mdPt,1:mdPt),'r')%Temporal
    plot(xSmpDegPos,rgc(mdPt:end,mdPt)','b')%Inferior
    plot(xSmpDegPos,rgc(mdPt,mdPt:end),'g')%Nasal
    plot(xSmpDegNeg,rgc(1:mdPt,mdPt)','k')%Superior 
    legend('Temporal','Inferior','Nasal','Superior')
    xlabel('Eccentricity (deg)'); set(gca,'XScale','log');
    ylabel('Denstiy (deg^{-2}');  set(gca,'YScale','log');

end 


end

