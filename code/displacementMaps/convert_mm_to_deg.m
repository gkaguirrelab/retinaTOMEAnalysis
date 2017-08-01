function ecc_deg = convert_mm_to_deg(ecc_mm)
%convert_mm_to_deg -- Converts milimeters to degrees on the retina based on
% the equation from the appendix of Watson 2014.
%
% Input: 
%   ecc_mm  = Retinal postion(s) in milimeters. Either scalar value or
%                     vector input accepted.
%
% Output: 
%   ecc_deg = Retinal postion(s) in degrees
%
% Sample Call:
%   ecc_mm         = 0:2:20;
%   sampleBaseDeg = convert_deg_to_mm(ecc_mm);
%
% MAB 2016

ecc_deg = 3.556.*ecc_mm + 0.05593.*(ecc_mm.^2) - 0.007358.*(ecc_mm.^3) +0.0003027.*(ecc_mm.^4);

end

