function [ data_mmSq ] = convert_degSq_to_mmSq( sampleBaseDeg, data_degSq )


% conversion from mm^2 to deg^2 (mm^2/deg^2)
% equation from end of Watson 2014 fit from Drasdo 1974

% Calculate the alpha conversion factor. It varies by eccentricity position
alpha= 0.0752+5.846e-5*sampleBaseDeg-1.064e-5*sampleBaseDeg.^2+4.116e-8*sampleBaseDeg.^3;

data_mmSq = data_degSq./alpha';

end

