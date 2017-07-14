function data_degSq = convert_mmSq_to_degSq(sampleBaseDeg, data_mmSq )
%convert_mmSq_to_degSq -- Converts square mm to square degrees
%
% data_mmSq = convert_mmSq_to_degSq( sampleBaseDeg, data_degSq )
%
% Usage: 
%   conversion from mm^2 to deg^2. 
%   equation from end of Watson 2014 fit from Drasdo 1974
% 
% Inputs:
%   sampleBaseDeg - sample position in degrees 
%   data_

% Calculate the alpha conversion factor. It varies by eccentricity position
alpha= 0.0752+5.846e-5*sampleBaseDeg-1.064e-5*sampleBaseDeg.^2+4.116e-8*sampleBaseDeg.^3;

data_degSq = data_mmSq.*alpha;

end

