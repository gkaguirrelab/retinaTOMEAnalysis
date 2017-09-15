function [ midgetRFDensity_degSq ] = getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, angle)
% getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, angle)
%
% This function implements Equation 8 of Watson 2014.
%
% Inputs
%   supportPosDeg - the positions (in degrees of visual angle) from the
%       fovea at which to calculate the midget receptive field density
%   angle - the meridian to evaluate. Acceptable values are: 
%       (0=nasal;90=superior;180=temporal;270=inferior)
%
% Outputs
%   midgetRFDensity_degSq - the density (receptive fields per square
%       degree) of midget receptive fields at each of the positions

%% Check the input

if sum([0 90 180 270]==angle) ~= 1
    error('The Watson equation for mRF density is defined only for the cardinal meridia');
end

%% Obtain the parameters of the modeled receptive field denisty for this angle

% Taken from Table 1 of Watson
switch angle
    case 0 % nasal
        a = 0.9729;
        r2 = 1.084;
        re = 7.633;        
        
    case 90 % superior
        a = 0.9935;
        r2 = 1.035;
        re = 16.35;        
        
    case 180 % temporal
        a = 0.9851;
        r2 = 1.058;
        re = 22.14;        
        
    case 270 % inferior
        a = 0.996;
        r2 = 0.9932;
        re = 12.13;        
        
end
        
% Cone density (cones / deg^2) at the foveal peak, same for every meridian
dcZero = 14804.6;

% Eccentricity (in degrees) at which the fraction of midget retinal
% ganglion cells has dropped by half from the initial fraction. The midget
% fraction as a function of eccentricity is assumed to be the same for all
% meridia
rm = 41.03;

% Note, the 2 * dcZero value is an expression for the number of midget
% retinal ganglion cells at the fovea, which is assumed to be 2 * dcZero.
% That is, it is assumed that the number of midget RGCs at the fovea is
% exactly equal to twice the number of cones.

midgetRFDensity_degSq = 2 * dcZero .* ...                                           % The number of RFs at the fovea (which is twice the cone denisty)
    ( (a*((1+supportPosDeg./r2)).^-2) + (1-a).*exp(-1.*supportPosDeg./re) ) .* ...  % The number of receptive fields as a function of eccentricity
    ((1+supportPosDeg./rm).^-1);                                                    % The midget fraction at each eccentricity


end

