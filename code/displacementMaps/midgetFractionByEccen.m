function midgetFraction = midgetFractionByEccen(ecc_deg)
%midgetFractionByEccen -- This function returns the fraction of midget 
% retinal ganglion cells in relation to the entire count of RGCs at a given
% distance (in deg) from the fovea. The equation is taken from Watson JoV 2014.   
%
% Inputs: 
%   ecc_deg        = Sample positions along a meridian in degrees.
%
% Outputs:
%   midgetFraction = Fraction of RGCs that are midget RGCs per retinal location. 
%
% Sample Call:
%   midgetFraction = midgetFractionByEccen(ecc_deg)
%
% MAB 2017

% Constants from Watson 2014.
f0 = 0.8928;
rm = 41.03;
% Equation for mRGC fraction.
midgetFraction = f0.*(1+(ecc_deg./rm)).^-1;
end 
