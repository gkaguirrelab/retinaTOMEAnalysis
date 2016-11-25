function [meridian]= curcio_4meridian(radDeg,smpPerDeg)

%this processes the curcio RGCd 4 meridian data

load('curcio_4meridian.mat')

ecc_mm = data(:,1);
ecc_deg = 3.556*ecc_mm+0.05992*ecc_mm.^2-0.007358*ecc_mm.^3+0.0003027*ecc_mm.^4;
alpha= 0.0752+5.846e-5*ecc_mm-1.064e-5*ecc_mm.^2+4.116e-8*ecc_mm.^3;
temp_deg2 = data(:,2).*(1./alpha);
sup_deg2 = data(:,4).*(1./alpha);
nasal_deg2 = data(:,6).*(1./alpha);
inferior_deg2 = data(:,8).*(1./alpha);

% Conversion of eccentricities in millimeters to degrees; Watson 2014

%%

[curve_nasal] = fit(ecc_deg,nasal_deg2 ,'smoothingspline','Exclude', find(isnan(nasal_deg2)),'SmoothingParam', 1);

[curve_sup] = fit(ecc_deg,sup_deg2,'smoothingspline', 'Exclude',find(isnan(sup_deg2)),'SmoothingParam', 1);

[curve_temp] = fit(ecc_deg,temp_deg2,'smoothingspline', 'Exclude',find(isnan(temp_deg2)),'SmoothingParam', 1);

[curve_inferior] = fit(ecc_deg,inferior_deg2,'smoothingspline', 'Exclude',find(isnan(inferior_deg2)),'SmoothingParam', 1);

meridian = zeros(2*(radDeg*smpPerDeg)+1,2*(radDeg*smpPerDeg)+1);
temp_hm = curve_temp(0:1/smpPerDeg:radDeg);
sup_vm = curve_sup(radDeg:-1/smpPerDeg:0)';
nasal_hm = curve_nasal(radDeg:-1/smpPerDeg:0);
inf_vm = curve_inferior(0:1/smpPerDeg:radDeg)';

[line,meridian(:,:,2),meridian(:,:,3)] = createGrid(radDeg,smpPerDeg);


%%interpolate 
for i = 1:size(meridian,1);
    for ii = 1:size(meridian,2)
        
        if meridian(i,ii,2) >= 0 & meridian(i,ii,2) <= 90
            VMd = curve_sup(meridian(i,ii,3));
            HMd = curve_nasal(meridian(i,ii,3));
            theta1 = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([0,90],[HMd,VMd],theta1);
            
        elseif meridian(i,ii,2) > 90 & meridian(i,ii,2) <= 180
            VMd = curve_sup(meridian(i,ii,3));
            HMd = curve_temp(meridian(i,ii,3));
            theta1 = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([90,180],[VMd,HMd],theta1);
            
        elseif meridian(i,ii,2) >= 180 & meridian(i,ii,2) <= 270
            VMd = curve_inferior(meridian(i,ii,3));
            HMd = curve_temp(meridian(i,ii,3));
            theta1 = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([180,270],[HMd,VMd],theta1);
            
        elseif meridian(i,ii,2) >= 270 & meridian(i,ii,2) <=360
            VMd = curve_inferior(meridian(i,ii,3));
            HMd = curve_nasal(meridian(i,ii,3));
            theta1 = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([270,360],[VMd,HMd],theta1);
            
        end
    end
end


mask=(meridian(:,:,3)<=radDeg);
double(mask(mask == 0)) = nan;
meridian(:,:,1) =meridian(:,:,1).*mask;

end

