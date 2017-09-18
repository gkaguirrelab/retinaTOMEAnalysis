function [ coneDensityPerSqDeg ] = calcTylerConeDensityByEccen( supportPosDeg, varargin )
% calcTylerConeDensityByEccen( supportPosDeg, varargin )
%
% This function implements the model for cone denisty proposed in:
%
%	Tyler, Christopher W. "Analysis of human receptor density."
%       Basic and clinical applications of vision science.
%       Springer Netherlands, 1997. 63-71.
%
% The motivating theory is that cone density (adjusted for cone surface
% area) is proportional to the incident light at any point on the retina. A
% distinctive feature of this model is that it captures the marked increase
% in cone density that occurs extreme eccentricities in proximity to the
% ora serrata.
%
% Retinal illuminance decreases with distance from the optical axis of the
% eye as incident light encounters the pupil aperture at an angle, and thus
% reduces the circular pupil area to an ellipse, eventually reaching a thin
% slit at the extreme periphery. This effect is partially compensated for
% by the decreased optical path length between the pupil and the retina at
% increasing eccentricities. Models of these two factors (in a model eye)
% provide for the receptorIlluminance as a function of eccentricity.
%
% Cone diameters increase as a function of eccentricity. This is modeled,
% and combined with the receptorIlluminance to provide luminousFlux.
%
% Cone density is finally modeled as the inverse of luminousFlux, subject
% to a scaling factors that sets peak cone density at the fovea.
%
% Implementation notes:
% - Equation (4) expresses the constant "u" as pi/2. This might lead to
%   the mistaken assumption that this constant is serving a trigonometric
%   role. In fact, the expression is in terms of degrees, and the constant
%   is just what was needed to provide an adequate fit of an arbitrary
%   function to the plots of Drasdo & Fowler 1974. This constant (which
%   I call here "anglePathLengthConstant") is set to 1.52, to emphasize
%   that it is not a trigonometric value.
% - The value of eyePathLengthConstant (w in the paper) is given as 0.19.
%   but I can only replicate figure 3B using a value of 0.125, found by
%   trial and error.
% - I have assumed that all of the original equations are performing
%   trigonometric operations in degrees.
% - Tyler identifies the innerSegDiamEccConstantDeg as the parameter that
%   is available to be adjusted to fit particular data sets (e.g.,
%   variation in density by meridian).
% - The equations are in density units of counts /  mm2. The routine takes
%   the eccentricity support in degrees, converts to mm internally, and
%   then converts back to degrees when returning the output
% - Interestingly, the rise in cone density at far eccentricies that is
%   apparent in units of mm^2 is attenuated in a plot of density in units
%   of deg^2.
%
% INPUTS (required)
%
% supportPosDeg - support for the returned function in units of degrees of
%   eccentricity on the retina.
%
% INPUTS (optional)
%   peakDensityDegSq - cone density at the fovea in units of counts / deg2.
%   innerSegDiamEccConstantDeg - a constant that controls the point at
%       which inner segment diameter starts to rise. Referred to as the
%       foveolar singularity constant.
%   eyePathLengthConstant - a constant related to the length of the model
%       eye that controls the optical path lenght model.
%   coneDiameterBase - controls the base size (related to the minimum) cone
%      diameter in the model, in units of micrometers.
%   coneDiameterExponent - the exponent of the cone diameter model
%   pupilAreaConstant - value related to the change in the effective area
%      of the pupil as a function of eccentricity.
%
% OUTPUTS
%
% coneDensityPerSqDeg - cone density in counts / deg2 as a function of
%   eccentricity in degrees, of the same length as the passed supportPosDeg
%   vector.
%
% Example call
%
%   close all
%   supportPosDeg = 0:0.01:90;
%   coneDensityPerSqDeg = calcTylerConeDensityByEccen( supportPosDeg, 'diagnosticFigures', true );
%   figure
%   loglog(supportPosDeg, coneDensityPerSqDeg);
%   xlabel('eccentricity [deg]');
%   ylabel('Cone density [count / deg2]');
%
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('supportPosDeg',@isnumeric);

% Optional
p.addParameter('diagnosticFigures',false,@islogical);
p.addParameter('peakDensityDegSq',14804.6,@isnumeric);
p.addParameter('innerSegDiamEccConstantDeg',0.2,@isnumeric);
p.addParameter('eyePathLengthConstant',0.125,@isnumeric);
p.addParameter('coneDiameterBase',1.25,@isnumeric);
p.addParameter('coneDiameterExponent',0.33,@isnumeric);
p.addParameter('pupilAreaConstant',0.82,@isnumeric);
p.addParameter('anglePathLengthConstant',1.52,@isnumeric);

%% Parse the parameters
p.parse(supportPosDeg, varargin{:});

% The scaling parameter reflects the cone density at the fovea, but in
% practice does not exactly match this number
scalingParam = p.Results.peakDensityDegSq ./ convert_mmSq_to_degSq(0,1) ;
innerSegDiamEccConstantDeg = p.Results.innerSegDiamEccConstantDeg;
eyePathLengthConstant = p.Results.eyePathLengthConstant;
coneDiameterBase = p.Results.coneDiameterBase;
coneDiameterExponent = p.Results.coneDiameterExponent;
pupilAreaConstant = p.Results.pupilAreaConstant;
anglePathLengthConstant = p.Results.anglePathLengthConstant;

% Perform the calculation
pupilArea = cos(degtorad(pupilAreaConstant.* supportPosDeg));
pathLength = 1 - (eyePathLengthConstant.*(1-cos(degtorad(anglePathLengthConstant.*supportPosDeg))));
receptorIlluminance = pupilArea ./ (pathLength.^2);
coneDiameters = coneDiameterBase .* (innerSegDiamEccConstantDeg + supportPosDeg).^coneDiameterExponent;
luminousFlux = 2 * pi .* coneDiameters.^2 .* receptorIlluminance;
coneDensityPerSqMM = scalingParam ./ luminousFlux;

% Convert mm^2 back to deg^2
degSqPerMmSq = convert_mmSq_to_degSq(supportPosDeg, ones(1,length(supportPosDeg) ));
coneDensityPerSqDeg = coneDensityPerSqMM .* degSqPerMmSq;

% Make some diagnostic figures if requested
if p.Results.diagnosticFigures
    
    % Replicate Tyler's Figure 3, showing the effect of retinal
    % eccentricity on the apparent area of the pupil and the reduction in
    % area of a the retinal image of a square degree of the visual world
    figure
    subplot(3,1,1);
    plot(supportPosDeg,pupilArea);
    ylim([0 1]);
    xlim([min(supportPosDeg) max(supportPosDeg)]);
    xlabel('eccentricity [deg]');
    ylabel('Relative pupil area');
    title('Figure 3A of Tyler 1997');
    subplot(3,1,2);
    plot(supportPosDeg,pathLength);
    ylim([0 1]);
    xlim([min(supportPosDeg) max(supportPosDeg)]);
    xlabel('eccentricity [deg]');
    ylabel('Relative retinal projection distance');
    title('Figure 3B of Tyler 1997');
    subplot(3,1,3);
    plot(supportPosDeg,receptorIlluminance);
    ylim([0 1]);
    xlim([min(supportPosDeg) max(supportPosDeg)]);
    xlabel('eccentricity [deg]');
    ylabel('Relative receptor illuminance');
    title('Figure 3B (inset) of Tyler 1997');
    hold off
    
    % Plot cone density in mm^2 as a function of eccentricity
    figure
    loglog(supportPosDeg, coneDensityPerSqMM / 1000,'.-');
    hold off
    xlabel('eccentricity [deg]');
    ylabel('Cone density (count / mm2) thousands');
    title('Figure 5 of Tyler 1997');
    
end % Make plots

end