% Examine the point of equivalence in the cumulative sums of the empirical
% RGC density and RF density

% Housekeeping

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
ecc_mm  = data(:,1);
ecc_deg = convert_mm_to_deg(ecc_mm);

rgcDensity_degSq_temporal = convert_mmSq_to_degSq(ecc_deg,data(:,2));
rgcDensity_degSq_superior = convert_mmSq_to_degSq(ecc_deg,data(:,4));
rgcDensity_degSq_nasal    = convert_mmSq_to_degSq(ecc_deg,data(:,6));
rgcDensity_degSq_inferior = convert_mmSq_to_degSq(ecc_deg,data(:,8));


% Adjust the RGC counts to reflect the fraction that is midget cells
midgetFraction = midgetFractionByEccen(ecc_deg);
midget_rgcDensity_degSq_temporal = rgcDensity_degSq_temporal.* midgetFraction;
midget_rgcDensity_degSq_superior = rgcDensity_degSq_superior.* midgetFraction;
midget_rgcDensity_degSq_nasal = rgcDensity_degSq_nasal.* midgetFraction;
midget_rgcDensity_degSq_inferior = rgcDensity_degSq_inferior.* midgetFraction;

% Plot the raw RGC density data
plot(ecc_deg,midget_rgcDensity_degSq_superior,'.r')
hold on
plot(ecc_deg,midget_rgcDensity_degSq_superior,'-k')
xlabel('eccentricity [deg]');
ylabel('density [cells / deg2]');
hold off

% Fit a spline to the RGC density data
splineFunctionSuperior = fit(ecc_deg,midget_rgcDensity_degSq_superior,'smoothingspline', 'Exclude',find(isnan(midget_rgcDensity_degSq_superior)),'SmoothingParam', 1);

% Switch now to equal eccentricity sampling
ecc_deg = 0:1:20;


%% Now get the RF density (reflecting the midget component)

% Receptive field density (Watson 2014, eq XX)
superiorRFDensity       = 2*(14804.6) * ( 0.9935*(1+ecc_deg/(1.035)).^-2+(1-0.9935)*exp(-1*ecc_deg/16.35));

% Adjust to have just the midget fraction component of the RF density
% (Watosn 2014 eq YY)
superiorRFDensity       = (superiorRFDensity .* (0.8928*(1+ecc_deg./41.03).^-1));

% Plot the raw RF density calculation
figure
plot(ecc_deg,superiorRFDensity,'.r')
hold on
plot(ecc_deg,superiorRFDensity,'-k')
xlabel('eccentricity [deg]');
ylabel('density [RF / deg2]');
hold off

% Watson 2*pi*r correction
ringArea = [0,diff(ecc_deg.^2 * pi)];


% Plot the RGC counts per ring
figure
plot(ecc_deg,splineFunctionSuperior(ecc_deg).*ringArea','.r')
hold on
plot(ecc_deg,splineFunctionSuperior(ecc_deg).*ringArea','-k')
xlabel('eccentricity [deg]');
ylabel('RGC counts per ring');
hold off


countPerRingRF = cumsum(superiorRFDensity.*ringArea);
countPerRingRGC = cumsum(splineFunctionSuperior(ecc_deg).*ringArea');

% Plot the RGC and RF density data, within growing areas
figure
plot(ecc_deg,countPerRingRGC,'.r')
hold on
plot(ecc_deg,countPerRingRGC,'-k')
xlabel('eccentricity [deg]');
ylabel('counts [cells per sector]');
plot(ecc_deg,countPerRingRF,'.r')
hold on
plot(ecc_deg,countPerRingRF,'-b')
hold off




function midgetFraction = midgetFractionByEccen(ecc_deg)
% function midgetFraction = midgetFractionByEccen(ecc_deg)
%
% This function returns the fraction of midget retinal ganglion cells in
% relation to the entire count of RGCs at a given distance (in deg) from the
% fovea. The equation is taken from Watson JoV 2016, 

f0 = 0.8928;
rm = 41.03;
midgetFraction = f0.*(1+(ecc_deg/rm)).^-1;
end % function

