%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scirpt that is setting us up to compute the displacement. 
% PARTS: 
%   *Fit Functions  - fits a frechet PDF to the RGC and fits a 2 term
%                    exponential to the RF and returns the 
%
%   *Plot the fits - Examine the fits ## COULD BE TURNED INTO A VERBOSE FLAG FOR THE FIT FUNCTIONS ## 
%
%   *Compute the number of receptive fields per 1 deg expanding ring
%
%   *Compute the number of midget retinal ganglion cells 1 deg expanding ring
%
%   *Calculate Displacement - NEEDs TO BE DONE WITH THE NEW FUNCTIONS FOR
%                             COUNT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%% Fit Functions
% Set the angle of the merdidian for the displacement to be computed on (0 = nasal, 90= superior).
angle = 0;
% Fit the RGC density with the Frechet PDF
[supportPosDeg,outParams_RGC,RGCdensityFit, scaleData] = fitRGCdensity(angle);
% Fit the RF density with the two term exponential 
RFfit = fitRFdensity(supportPosDeg,angle,scaleData);

%% Plot the fits
figure;
plot(supportPosDeg,RGCdensityFit,'g')
hold on
plot(supportPosDeg,RFfit(supportPosDeg),'r')


%% compute the number of receptive fields per 1 deg expanding ring 
verbose = true;
[CountPerRingRF] = RFcountFunc(RFfit,verbose,scaleData)


%% compute the number of midget retinal ganglion cells per 1 deg expanding ring 
verbose = true;
[CountPerRingRGC] = RGCcountFunc(outParams_RGC,verbose,scaleData)

%% Calculate Displacement 
% NEEDS TO BE DONE FOR THE NEW METHOD OF COUNTS PER RING.
