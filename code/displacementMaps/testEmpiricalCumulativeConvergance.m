% #############################################################
% Script to "reproduce" (get close) Figure 1 from Turpin/Mckendrick 2015.
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

clc
close all

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
supportPosMm  = data(:,1);
supportPosDeg = convert_mm_to_deg(supportPosMm);

rgcDensity_degSq_temporal = convert_mmSq_to_degSq(supportPosDeg,data(:,2));
rgcDensity_degSq_superior = convert_mmSq_to_degSq(supportPosDeg,data(:,4));
rgcDensity_degSq_nasal    = convert_mmSq_to_degSq(supportPosDeg,data(:,6));
rgcDensity_degSq_inferior = convert_mmSq_to_degSq(supportPosDeg,data(:,8));


% Adjust the RGC counts to reflect the fraction that is midget cells
midgetFraction = midgetFractionByEccen(supportPosDeg);
midget_rgcDensity_degSq_temporal = rgcDensity_degSq_temporal.* midgetFraction;
midget_rgcDensity_degSq_superior = rgcDensity_degSq_superior.* midgetFraction;
midget_rgcDensity_degSq_nasal = rgcDensity_degSq_nasal.* midgetFraction;
midget_rgcDensity_degSq_inferior = rgcDensity_degSq_inferior.* midgetFraction;

% Plot the raw RGC density data
plot(supportPosDeg,midget_rgcDensity_degSq_superior,'.r')
hold on
plot(supportPosDeg,midget_rgcDensity_degSq_superior,'-k')
xlabel('eccentricity [deg]');
ylabel('density [cells / deg2]');
hold off

% Fit a spline to the RGC density data
splineFunctionSuperior = fit(supportPosDeg,midget_rgcDensity_degSq_superior,'smoothingspline', 'Exclude',find(isnan(midget_rgcDensity_degSq_superior)),'SmoothingParam', 1);

% Switch now to equal eccentricity sampling
supportPosDeg = 0:1:20;


%% Now get the RF density (reflecting the midget component)

% Receptive field density (Watson 2014, eq XX)
superiorRFDensity       = 2*(14804.6) * ( 0.9935*(1+supportPosDeg/(1.035)).^-2+(1-0.9935)*exp(-1*supportPosDeg/16.35));

% Adjust to have just the midget fraction component of the RF density
% (Watosn 2014 eq YY)
superiorRFDensity       = (superiorRFDensity .* (0.8928*(1+supportPosDeg./41.03).^-1));

% Plot the raw RF density calculation
figure
plot(supportPosDeg,superiorRFDensity,'.r')
hold on
plot(supportPosDeg,superiorRFDensity,'-k')
xlabel('eccentricity [deg]');
ylabel('density [RF / deg2]');
hold off

% Watson 2*pi*r correction I AM NOT SURE 
ringArea = [0,diff(supportPosDeg.^2 * pi)];


% Plot the RGC counts per ring
figure
plot(supportPosDeg,splineFunctionSuperior(supportPosDeg).*ringArea','.r')
hold on
plot(supportPosDeg,splineFunctionSuperior(supportPosDeg).*ringArea','-k')
xlabel('eccentricity [deg]');
ylabel('RGC counts per ring');
hold off


countPerRingRF = cumsum(superiorRFDensity.*ringArea);
countPerRingRGC = cumsum(splineFunctionSuperior(supportPosDeg).*ringArea');

% Plot the RGC and RF density data, within growing areas
figure
plot(supportPosDeg,countPerRingRGC,'.r')
hold on
plot(supportPosDeg,countPerRingRGC,'-k')
xlabel('eccentricity [deg]');
ylabel('counts [cells per sector]');
plot(supportPosDeg,countPerRingRF,'.r')
hold on
plot(supportPosDeg,countPerRingRF,'-b')
hold off




