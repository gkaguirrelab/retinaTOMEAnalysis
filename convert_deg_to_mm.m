function [ sampleBase_mm ] = convert_deg_to_mm( sampleBase_deg )

% the eccentricy of each sample in MM
% equation from end of Watson 2014
sampleBase_mm = 0.268.*sampleBase_deg + 0.0003427.*(sampleBase_deg.^2) - 8.3309e-6.*(sampleBase_deg.^3);

end

