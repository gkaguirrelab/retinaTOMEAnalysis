function CountPerRingRGC = RGCcountFunc(outParams_RGC,verbose)
%RGCcountFunc -- returns a function that calculates the number of RGCs in a 1 deg ring  
% 
% Decription:
%   This returns a function that computes the definite intergral of the RGC density function * 
%   2*pi*r from r to r+1 deg. This give the number of cell in that 1 deg
%   ring statring a r and ending at r+1. 
%
% Inputs:
%   outParams_RGC   = The params fit from fitRGCdensity.m. These are the
%                     parameters of the frechet fit 
%   verbose         = Option to plot the output function.
%
% Outputs:
%   CountPerRingRGC = Function that estimates the total number of RGCs in a 1 deg ring
%
% Sample Call: 
%   angle   = 0;
%   verbose = true;
%   [supportPosDeg,outParams_RGC,RGCdensityFit, scaleData] = fitRGCdensity(angle);
%   CountPerRingRGC = RGCcountFunc(outParams_RGC,verbose)

% MAB 2017

% Parameters fit from the frechet fit to the curcio 1990 data of RCG density 
a    = outParams_RGC(1);
s    = outParams_RGC(2);
m    = outParams_RGC(3);

% The definite integral of 2*pi*r * Frechet PDF. This should give the
% number of cells in an expanding 1 deg ring. 
CountPerRingRGC = @(x)(2.*pi.*((x.*exp(-(((x+1)-m)/s).^-a) - (a.^-1).*(((x+1)-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,(((x+1)-m)/s).^a)))) - ...
    2.*pi.*((x.*exp(-((x-m)/s).^-a) - (a.^-1).*((x-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,((x-m)/s).^a)))));

% Plot Option 
if verbose == true
    figure
    subplot(1,2,1)
    plot(1:20,CountPerRingRGC(1:20).*scaleData)
    title('RGC Cell Count')
    xlabel('eccentricity (deg)')
    ylabel('RGC Count (cells)')
    subplot(1,2,2)
    plot(1:20,cumsum(CountPerRingRGC(1:20).*scaleData))
    title('Cumulative RGC Cell Count')
    xlabel('eccentricity (deg)')
    ylabel('Cumulative RGC Count (cells)')
end