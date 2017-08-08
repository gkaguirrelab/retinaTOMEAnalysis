function CountPerRingRF = RFcountFunc(RFfit,verbose,scaleData)
%RFcountFunc -- returns a function that calculates the number of RFs in a 1 deg ring  
% 
% Decription:
%   This returns a function that computes the definite intergral of the RF density function * 
%   2*pi*r from r to r+1 deg. This give the number of RFs in that 1 deg
%   ring statring at r and ending at r+1. 
%
% Inputs:
%   outParams_RGC   = The params fit from fitRFdensity.m. These are the
%                     parameters of the two term exponential. 
%   verbose         = Option to plot the output function.
% Outputs:
%   CountPerRingRF = Function that estimates the total number of RFs in a 1 deg ring
%
% Sample Call: 
%   RFfit = fitRFdensity(supportPosDeg,angle,scaleData);
%   [CountPerRingRF] = RFcountFunc(RFfit,verbose)

% MAB 2017


% Parameters fit to a two term expontial of the form f(x) = a*exp(b*x) + c*exp(d*x) 
% b and d are *-1 becuase of a mismatch bewteen the way matlab writes the
% equation and the way i did the math for the integral but this should not
% be a problem...

a = RFfit.a;
b = RFfit.b.*-1;
c = RFfit.c;
d = RFfit.d.*-1;

% the definite integral of 2*pi*r * two-term exponential decay function. This should give the
% number of recptive fields in an expanding 1 deg ring. 
CountPerRingRF = @(x)(2.*pi.*((x+1).*(-a.*exp(-b.*(x+1))/b - c.*exp(-d.*(x+1))/d) - (exp(-b.*(x+1))/(b.^2) + c.*exp(-d.*(x+1))/(d.^2))) - ...
    2.*pi.*(x.*(-a.*exp(-b.*x)/b - c.*exp(-d.*x)/d) - (exp(-b.*x)/(b.^2) + c.*exp(-d.*x)/(d.^2))));


%Plot Option
if verbose == true
    figure
    subplot(1,2,1)
    plot(1:0.2:20,CountPerRingRF(1:0.2:20).*scaleData)
    title('RF Count')
    xlabel('eccentricity (deg)')
    ylabel('RF Count (receptive fields)')
    subplot(1,2,2)
    plot(1:0.2:20,cumsum(CountPerRingRF(1:0.2:20).*scaleData))
    title('Cumulative RF Count')
    xlabel('eccentricity (deg)')
    ylabel('Cumulative RF Count (receptive fields)')
end