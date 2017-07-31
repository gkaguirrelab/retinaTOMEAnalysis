function sampleBaseDeg = convert_mm_to_deg(sampleBaseMm)
%convert_mm_to_deg -- Converts milimeters to degrees on the retina based on
% the equation from the appendix of Watson 2014.
%
% Input: 
%   sampleBaseMm   = Retinal postion(s) in milimeters. Either scalar value or
%                     vector input accepted.
%
% Output: 
%   sampleBaseDeg  = Retinal postion(s) in degrees
%
% Sample Call:
%   sampleBaseMm   = 0:2:20;
%   sampleBaseMDeg = convert_deg_to_mm(sampleBaseDeg);
%
% MAB 2016

sampleBaseDeg = 3.556.*sampleBaseMm + 0.05593.*(sampleBaseMm.^2) - 0.007358.*(sampleBaseMm.^3) +0.0003027.*(sampleBaseMm.^4);

end

