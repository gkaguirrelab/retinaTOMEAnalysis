% #############################################################
% Script to reproduce Figure 1 from Turpin/Mckendrick 2015.
% 
% The goal is to produce a plot of the cumulative retinal
% ganglion cell and receptive field counts as a function of
% eccentricity (mm). 
% 
% Everthing should be in units of mm, mm^2, or cell count
% Conversion between units found at the end of Watson 2014
% and fit from Drasdo 1974.
% 
% MAB 2016
%
% Notes:
% Curcio and Allen 1990 reference frames is the retinal feild
%       Data loaded from .xml is in the retinal coordinate frame
%       No extra steps neded to reference the data.
% Drasdo 2007 usues the visual field refrerence frame
%       Indexing the output of densityRF():
%       ex. RFdensity = densityRf(radDeg,smpPerDeg,'OFF');
%            midPoint = center pixel of output image
%            Temporal retina = left of midPoint
%            Nasal Retina    = right of midPoint
%            Supior Retina   = above midPoint
%            Inferior Retina = below midPoint
% Turpin McKendrick 2015 reference frame is the visual field
%
% Watson 2014 reference frame is the visual field
%
% THIS CODE SETS THE REFERENCE FRAME TO THE RETINAL FIELD
%
% #############################################################

%% Load the RGC Density Data from Curcio and Allen 1990:
% The meridian assignments are in the retinal coordinate frame. This can be
% verified by observing that there is an interruption in the count data for
% the nasal meridian corresponding to the blind spot.
load('curcio_4meridian.mat') % Load the data 
sampleBase_RGC_mm  = data(:,1); % Assign the eccentrciy (mm) to a var
RGCdenisty_mmSq_temporal = data(:,2); % Assign the temporal RGC denstiy (cells/mm^2) to a var
RGCdenisty_mmSq_superior  = data(:,4); % Assign the superior RGC denstiy (cells/mm^2) to a var
RGCdenisty_mmSq_nasal  = data(:,6); % Assign the nasal RGC denstiy (cells/mm^2) to a var
RGCdenisty_mmSq_inferior  = data(:,8); % Assign the inferior RGC denstiy (cells/mm^2) to a var

%% Generate the Receptive Field Density From the Drasdo 2007:
% Set Parameters 
radDeg      = 20; % Radius in Degrees
smpPerDeg   = 2; % Samples per Degree

% Obtain the Receptive field density per square degree within the sampling
% area specified (also in degrees of visual angle)
RFdensity_sqDeg = densityRf(radDeg,smpPerDeg,'full'); % Generates a 2D Receptive Field Density plot 

% Get data from Superior Merdian as a first pass check to validate the
% pipeline of Turpin/McKendrick 

midPoint = round(size(RFdensity_sqDeg,1)/2); % find middle of RF density image

% Grab a vector that corresponds to the superior retinal meridian
sup_RFdensity_sqDeg = RFdensity_sqDeg(1:midPoint,midPoint); % extract the superior meridian from 2D Receptive Field Desity plot 
sup_RFdensity_sqDeg = flipud(sup_RFdensity_sqDeg); % flip column vector so 0 deg is at top

% Calculate a sample base
sampleBase_RF_deg = (0:1/smpPerDeg:radDeg)'; % the eccentricy of each sample of sup_RFdensity in degrees 

% Convert the RF counts (and sample base) from degress / deg^2 to mm and cells/mm^2
sup_RFdensity_mm = convert_degSq_to_mmSq(sampleBase_RF_deg, sup_RFdensity_sqDeg);
sampleBase_RF_mm=convert_deg_to_mm(sampleBase_RF_deg);

%% Fit splines to the data RF and RGC data
[sup_RFdensity_mm_fit] = fit(sampleBase_RF_mm,sup_RFdensity_mm,'smoothingspline','Exclude', find(isnan(sup_RFdensity_mm)),'SmoothingParam', 1);
[sup_RGCdensity_mm_fit] = fit(sampleBase_RGC_mm,RGCdenisty_mmSq_superior,'smoothingspline','Exclude', find(isnan(RGCdenisty_mmSq_superior)),'SmoothingParam', 1);


%% Calculate sector size to extract cell count from cells/mm^2
sectorAngle = 6; % Angle of the sector ### turpin code does not use sector angle or pi

annulus_radius_mm = 0.005; % parameter from Turpin code

% center of the segment in mm eccentricity
radii_mm = 0:annulus_radius_mm:5; % vector of radii to match Turpin code

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
