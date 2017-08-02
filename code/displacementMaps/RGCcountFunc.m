function [CountPerRingRGC] = RGCcountFunc(outParams_RGC,verbose)

%RGCcountFunc -- the definite intergral of the RGC density function * 2*pi*r
% from r to r+1 UNITS UNSURE FOR +1
% This script
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%

% Parameters fit from the frechet fit to the curcio 1990 data of RCG density 
a    = outParams_RGC(1);
s    = outParams_RGC(2);
m    = outParams_RGC(3);


% the definite integral of 2*pi*r * Frechet PDF. This should give the
% number of cells in an expanding 1mm ring. 
CountPerRingRGC = @(x)(2.*pi.*((x.*exp(-(((x+1)-m)/s).^-a) - (a.^-1).*(((x+1)-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,(((x+1)-m)/s).^a)))) - ...
    2.*pi.*((x.*exp(-((x-m)/s).^-a) - (a.^-1).*((x-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,((x-m)/s).^a)))));

% Plot
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