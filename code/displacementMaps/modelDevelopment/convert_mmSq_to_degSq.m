function samplesDegSq = convert_mmSq_to_degSq(supportPosDeg, samplesMmSq )
%convert_mmSq_to_degSq --  Converts square milimeters to square degrees on 
% the retina based on the equation from the appendix of Watson 2014 which 
% is a fit to a plot in Drasdo and Fowler 1974.
%
% Inputs:
%   supportPosDeg = sample position in degrees 
%   samplesMmSq   = Counts per square mm at locations corresponding to the sampleBaseDeg 
%
% OutPuts:
%   samplesDegSq  = Counts per square deg at locations corresponding to the sampleBaseDeg 
%
% Sample Call
%   sampleBaseDeg = [0, 5, 10, 15, 20]:
%   samplesMmSq   = [10000, 1000, 500, 100, 10];
%   samplesDegSq  = convert_mmSq_to_degSq(sampleBaseDeg, samplesMmSq)
% 
% MAB 2016

% Calculate the alpha conversion factor. It varies by eccentricity position
% Units of alpha = mm^2/deg^2
alpha= 0.0752+5.846e-5*supportPosDeg-1.064e-5*supportPosDeg.^2+4.116e-8*supportPosDeg.^3;

% Apply alpha to convert mm^2 to deg^2  (Cells/mm^2 * mm^2/deg^2 = Cells/deg^2)
samplesDegSq = samplesMmSq.*alpha;

end

