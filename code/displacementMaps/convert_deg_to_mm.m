function ecc_mm = convert_deg_to_mm(ecc_deg)
%convert_deg_to_mm -- Converts degrees to millimeters on the retina based 
% on the equation from the appendix of Watson 2014.
%
% Input: 
%   ecc_deg = Retinal postion(s) in degrees. Either scalar value or
%             vector input accepted.
%
% Output: 
%   ecc_mm  = Retinal postion(s) in millimeters
%
% Sample Call:
%   sampleBaseDeg = 0:2:20;
%   sampleBaseMm  = convert_deg_to_mm(sampleBaseDeg);
%
% MAB 2016

ecc_mm = 0.268.*ecc_deg + 0.0003427.*(ecc_deg.^2) - 8.3309e-6.*(ecc_deg.^3);

end

