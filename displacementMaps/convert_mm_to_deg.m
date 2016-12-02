function [ sampleBase_deg ] = convert_mm_to_deg( sampleBase_mm )

% the eccentricy of each sample in MM
% equation from end of Watson 2014
sampleBase_deg = 3.556.*sampleBase_mm + 0.05593.*(sampleBase_mm.^2) - 0.007358.*(sampleBase_mm.^3) +0.0003027.*(sampleBase_mm.^4);

end

