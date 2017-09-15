function [ coneDensityPerSqMM ] = tylerConeDensityPerSqMM( supportPosDeg, varargin )
% tylerConeDensityPerSqDeg( supportPosDeg, varargin )
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
% - It is difficult to determine from the printed equation (4) in
%   Tyler (1997) the operator that follows the constant "u". While I would
%   have presumed a "-", as this would correspond to taking the
%   complementary angle of theta, I can only replicate Figure 3B by making
%   this operator multiplication.
% - The value of eyePathLengthConstant (w in the paper) is given as 0.19.
%   but I can only replicate figure 3B using a value of 0.125, found by
%   trial and error.
% - I have assumed that all of the original equations are performing
%   trigonometric operations in degrees.
% - Tyler identifies the innerSegDiamEccConstantDeg as the parameter that
%   is available to be adjusted to fit particular data sets (e.g.,
%   variation in density by meridian).
% - The equations are in density units of counts /  mm2.
%
% INPUTS (required)
%
% supportPosDeg - support for the returned function in units of degrees of
%   eccentricity on the retina.
%
% INPUTS (optional)
%   peakDensityMMSq - cone density at the fovea in units of counts / mm2.
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
% coneDensityPerSqMM - cone density in counts / mm2 as a function of
%   eccentricity in degrees, of the same length as the passed supportPosDeg
%   vector.
%
% Example call
%
%   close all
%   supportPosDeg = 0:0.01:90;
%   coneDensityPerSqMM = tylerConeDensityPerSqMM( supportPosDeg, 'diagnosticFigures', true );
%   figure
%   loglog(supportPosDeg, coneDensityPerSqMM / 1000);
%   hold on
%   coneDensityPerSqMM = tylerConeDensityPerSqMM( supportPosDeg, 'peakDensityMMSq', 36500, 'innerSegDiamEccConstantDeg', 0.08 );
%   loglog(supportPosDeg, coneDensityPerSqMM / 1000,'.-');
%   hold off
%   xlabel('eccentricity [deg]');
%   ylabel('Cone density (count / mm2) thousands');
%   title('Figure 5 of Tyler 1997');
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('supportPosDeg',@isnumeric);

% Optional
p.addParameter('diagnosticFigures',false,@islogical);
p.addParameter('peakDensityMMSq',30000,@isnumeric);
p.addParameter('innerSegDiamEccConstantDeg',0.2,@isnumeric);
p.addParameter('eyePathLengthConstant',0.125,@isnumeric);
p.addParameter('coneDiameterBase',1.25,@isnumeric);
p.addParameter('coneDiameterExponent',0.33,@isnumeric);
p.addParameter('pupilAreaConstant',0.82,@isnumeric);

%% Parse the parameters
p.parse(supportPosDeg, varargin{:});

peakDensityMMSq = p.Results.peakDensityMMSq;
innerSegDiamEccConstantDeg = p.Results.innerSegDiamEccConstantDeg;
eyePathLengthConstant = p.Results.eyePathLengthConstant;
coneDiameterBase = p.Results.coneDiameterBase;
coneDiameterExponent = p.Results.coneDiameterExponent;
pupilAreaConstant = p.Results.pupilAreaConstant;
halfPi = pi/2;

% Perform the calculation
pupilArea = cos(degtorad(pupilAreaConstant.* supportPosDeg));
pathLength = 1 - (eyePathLengthConstant.*(1-cos(degtorad(halfPi.*supportPosDeg))));
receptorIlluminance = pupilArea ./ (pathLength.^2);
coneDiameters = coneDiameterBase .* (innerSegDiamEccConstantDeg + supportPosDeg).^coneDiameterExponent;
luminousFlux = 2 * pi .* coneDiameters.^2 .* receptorIlluminance;
coneDensityPerSqMM = (peakDensityMMSq * 2 * pi * 2 * pi ) ./ luminousFlux;

% Make some diagnostic figures if requested
if p.Results.diagnosticFigures
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
end

end