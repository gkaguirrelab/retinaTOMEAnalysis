function sampleBaseMm = convert_deg_to_mm(sampleBaseDeg)
%convert_deg_to_mm -- Converts degrees to millimeters on the retina based 
% on the equation from the appendix of Watson 2014.
%
% Input: 
%   sampleBaseDeg = Retinal postion(s) in degrees. Either scalar value or
%                   vector input accepted.
%
% Output: 
%   sampleBaseMm  = Retinal postion(s) in millimeters
%
% Sample Call:
%   sampleBaseDeg = 0:2:20;
%   sampleBaseMm  = convert_deg_to_mm(sampleBaseDeg);
%
% MAB 2016

sampleBaseMm = 0.268.*sampleBaseDeg + 0.0003427.*(sampleBaseDeg.^2) - 8.3309e-6.*(sampleBaseDeg.^3);

end

