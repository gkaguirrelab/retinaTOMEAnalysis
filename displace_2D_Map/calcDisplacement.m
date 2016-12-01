function Displacement = calcDisplacement(RFdensity_mm,sampleBase_RF_mm,RGCdenisty_mmSq,sampleBase_RGC_mm,radMM,smpPerMM,sectorAngle)
% Angle of the sector ### turpin code does not use sector angle or pi



%% Fit splines to the data RF and RGC data
[sup_RFdensity_mm_fit] = fit(sampleBase_RF_mm',RFdensity_mm','smoothingspline','Exclude', find(isnan(RFdensity_mm)),'SmoothingParam', 1);
[sup_RGCdensity_mm_fit] = fit(sampleBase_RGC_mm',RGCdenisty_mmSq','smoothingspline','Exclude', find(isnan(RGCdenisty_mmSq)),'SmoothingParam', 1);


%% Calculate sector size to extract cell count from cells/mm^2


annulus_radius_mm = 1/smpPerMM; % parameter from Turpin code

% center of the segment in mm eccentricity
radii_mm = 0:annulus_radius_mm:radMM; % vector of radii to match Turpin code

areaPerSeg_mmSq = ( (pi*(radii_mm + annulus_radius_mm/2).^2) - ...
                    (pi*(radii_mm - annulus_radius_mm/2).^2) ) * ...
                    (sectorAngle/360);

%% Apply the sectors to the fit Densities 
countRF = sup_RFdensity_mm_fit(radii_mm).*areaPerSeg_mmSq'; % multiply RFs/mm^2 *mm^2 to get RF count per sector 
countRGC = sup_RGCdensity_mm_fit(radii_mm).*areaPerSeg_mmSq';  % multiply RGCs/mm^2 *mm^2 to get RGC count per sector 

countRFsum = cumsum(countRF); % Cumulative sum of the recetive fields 
countRGCsum = cumsum(countRGC); % Cumulative sum of the recetive fields 

%% Make plot from Turpin/McKencdrick fig. 1
figure
hold on 
plot(radii_mm,countRFsum,'r')
plot(radii_mm,countRGCsum,'b')
legend('Receptive Fields','Retinal Ganglion Cell')
xlabel('Eccentricity (mm)')
ylabel('Cumulative RGC/RF Count')

% Calculate the displacement by finding the mm difference at equivalent
% count points
mmPerRGCcountAtRFcountPositions=interp1(countRFsum,radii_mm,countRGCsum,'spline');
Displacement=abs(radii_mm'-mmPerRGCcountAtRFcountPositions);

weights = radii_mm;
weights(weights <= 1.5) = 1;
weights(weights >1.5 ) = 0;

fitParams = fitGammaToDisplacement(radii_mm, Displacement', weights)

end