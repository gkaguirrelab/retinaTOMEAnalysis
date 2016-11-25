function [rgcDensity]= rgcDensityMM(radMM,smpPerMM)
%This function  the curcio RGCd 4 meridian data

%% Load and sort the data
load('curcio_4meridian.mat')
ecc_mm = data(:,1);
temp_mm = data(:,2);
sup_mm = data(:,4);
nasal_mm = data(:,6);
inferior_mm = data(:,8);

%% Fit the data with splines 
[curve_nasal] = fit(ecc_mm,nasal_mm ,'smoothingspline','Exclude', find(isnan(nasal_mm)),'SmoothingParam', 1);

[curve_sup] = fit(ecc_mm,sup_mm,'smoothingspline', 'Exclude',find(isnan(sup_mm)),'SmoothingParam', 1);

[curve_temp] = fit(ecc_mm,temp_mm,'smoothingspline', 'Exclude',find(isnan(temp_mm)),'SmoothingParam', 1);

[curve_inferior] = fit(ecc_mm,inferior_mm,'smoothingspline', 'Exclude',find(isnan(inferior_mm)),'SmoothingParam', 1);

%% generate empty output matrix 
rgcDensity(:,:,3) = zeros(2*(radMM*smpPerMM)+1,2*(radMM*smpPerMM)+1);

%% Create ecc and pol angle map to size
[~,rgcDensity(:,:,1),rgcDensity(:,:,2)] = createGrid(radMM,smpPerMM);


%% Interpolate points between the meridians  
for i = 1:size(rgcDensity,1);
    for ii = 1:size(rgcDensity,2)
        
        if rgcDensity(i,ii,1) >= 0 & rgcDensity(i,ii,1) <= 90
            VMd = curve_sup(rgcDensity(i,ii,2));
            HMd = curve_nasal(rgcDensity(i,ii,2));
            theta1 = rgcDensity(i,ii,1);
            rgcDensity(i,ii,3) = interp1([0,90],[HMd,VMd],theta1);
            
        elseif rgcDensity(i,ii,1) > 90 & rgcDensity(i,ii,1) <= 180
            VMd = curve_sup(rgcDensity(i,ii,2));
            HMd = curve_temp(rgcDensity(i,ii,2));
            theta1 = rgcDensity(i,ii,1);
            rgcDensity(i,ii,3) = interp1([90,180],[VMd,HMd],theta1);
            
        elseif rgcDensity(i,ii,1) >= 180 & rgcDensity(i,ii,1) <= 270
            VMd = curve_inferior(rgcDensity(i,ii,2));
            HMd = curve_temp(rgcDensity(i,ii,2));
            theta1 = rgcDensity(i,ii,1);
            rgcDensity(i,ii,3) = interp1([180,270],[HMd,VMd],theta1);  
            
        elseif rgcDensity(i,ii,1) >= 270 & rgcDensity(i,ii,1) <=360
            VMd = curve_inferior(rgcDensity(i,ii,2));
            HMd = curve_nasal(rgcDensity(i,ii,2));
            theta1 = rgcDensity(i,ii,1);
            rgcDensity(i,ii,3) = interp1([270,360],[VMd,HMd],theta1);      
        end
    end
end

%% Mask the output
mask=(rgcDensity(:,:,2)<=radMM);
double(mask(mask == 0)) = nan;
rgcDensity(:,:,3) =rgcDensity(:,:,3).*mask;

end

