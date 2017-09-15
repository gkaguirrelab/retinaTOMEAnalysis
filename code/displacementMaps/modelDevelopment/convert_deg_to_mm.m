function supportPosMm = convert_deg_to_mm(supportPosDeg)
%convert_deg_to_mm -- Converts degrees to millimeters on the retina based 
% on the equation from the appendix of Watson 2014.
%
% Input: 
%   supportPosDeg = Retinal postion(s) in degrees. Either scalar value or
%                   vector input accepted.
%
% Output: 
%   supportPosMm  = Retinal postion(s) in millimeters
%
% Sample Call:
%   supportPosDeg = 0:2:20;
%   supportPosMm  = convert_deg_to_mm(supportPosDeg);
%
% MAB 2016

supportPosMm = 0.268.*supportPosDeg + 0.0003427.*(supportPosDeg.^2) - 8.3309e-6.*(supportPosDeg.^3);

end

