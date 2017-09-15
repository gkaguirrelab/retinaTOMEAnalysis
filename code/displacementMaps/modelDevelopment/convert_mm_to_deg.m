function supportPosDeg = convert_mm_to_deg(supportPosMm)
%convert_mm_to_deg -- Converts milimeters to degrees on the retina based on
% the equation from the appendix of Watson 2014.
%
% Input: 
%   supportPosMm  = Retinal postion(s) in milimeters. Either scalar value or
%                     vector input accepted.
%
% Output: 
%   supportPosDeg = Retinal postion(s) in degrees
%
% Sample Call:
%   supportPosMm  = 0:2:20;
%   supportPosDeg = convert_deg_to_mm(supportPosMm);
%
% MAB 2016

supportPosDeg = 3.556.*supportPosMm + 0.05593.*(supportPosMm.^2) - 0.007358.*(supportPosMm.^3) +0.0003027.*(supportPosMm.^4);

end

