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
% The meridian assignments are in the retinal coordinate frame. This can be
% verified by observing that there is an interruption in the count data for
% the nasal meridian corresponding to the blind spot.
load('source_data/curcio_4meridian.mat')

ecc_mm = data(:,1);
temp_mmSq = data(:,2);
sup_mmSq = data(:,4);
nasal_mmSq = data(:,6);
inferior_mmSq = data(:,8);

% Conversion of eccentricities in millimeters to degrees; Watson 2014

%%

[curve_nasal] = fit(ecc_mm,nasal_mmSq,'smoothingspline','Exclude', find(isnan(nasal_mmSq)),'SmoothingParam', 1);

[curve_sup] = fit(ecc_mm,sup_mmSq,'smoothingspline', 'Exclude',find(isnan(sup_mmSq)),'SmoothingParam', 1);

[curve_temp] = fit(ecc_mm,temp_mmSq,'smoothingspline', 'Exclude',find(isnan(temp_mmSq)),'SmoothingParam', 1);

[curve_inferior] = fit(ecc_mm,inferior_mmSq,'smoothingspline', 'Exclude',find(isnan(inferior_mmSq)),'SmoothingParam', 1);

meridian = zeros(2*(round(radMM)*round(smpPerMM))+1,2*(round(radMM)*round(smpPerMM))+1);

[~,meridian(:,:,2),meridian(:,:,3)] = createGrid(radMM,smpPerMM);


%%interpolate 
for i = 1:size(meridian,1);
    for ii = 1:size(meridian,2)
        
        if meridian(i,ii,2) >= 0 & meridian(i,ii,2) <= 90
            VMd = curve_sup(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_nasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([0,90],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) > 90 & meridian(i,ii,2) <= 180
            VMd = curve_sup(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_temp(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([90,180],[VMd,HMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 180 & meridian(i,ii,2) <= 270
            VMd = curve_inferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_temp(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([180,270],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 270 & meridian(i,ii,2) <=360
            VMd = curve_inferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_nasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([270,360],[VMd,HMd],theta,interp);
            
        end
    end
end 

mask=(meridian(:,:,3)<=radMM);
double(mask(mask == 0)) = nan;
RGCdensity =meridian(:,:,1).*mask;
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

