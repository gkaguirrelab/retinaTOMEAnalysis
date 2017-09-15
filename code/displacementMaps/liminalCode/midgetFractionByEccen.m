function midgetFraction = midgetFractionByEccen(supportPosDeg)
% midgetFractionByEccen -- This function returns the fraction of midget 
% retinal ganglion cells in relation to the entire count of RGCs at a given
% distance (in deg) from the fovea.
%
% The equation is taken from Watson JoV 2014 (eq 7, plotted in fig 8), 
% which Watson took from Drasdo et al 2007.
%
% Inputs: 
%   supportPosDeg  = Sample positions along a meridian in degrees.
%
% Outputs:
%   midgetFraction = Fraction of RGCs that are midget RGCs per retinal location. 
%
% Sample Call:
%   midgetFraction = midgetFractionByEccen(supportPosDeg)
%
% MAB 2017

% Constants from Watson 2014.
f0 = 0.8928;
rm = 41.03;

% Equation for mRGC fraction.
midgetFraction = f0.*(1+(supportPosDeg./rm)).^-1;
end 
