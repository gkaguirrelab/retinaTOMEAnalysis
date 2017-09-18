function [rgcDensityFit] = getSplineFitToRGCDensity(polarAngle, varargin)
% getSplineFitToRGCDensity(angle)
%
% This routine returns a fit to RGC density data. Raw RGC density data
% are taken from Curcio et al (1990) unless otherwise specified.
% Fits to the cardinal meridians are obtained. Interpolation over the
% parameters of the fit are used to produce a function that returns RGC
% density for an arbitrary meridian angle.
%
% Inputs:
%   polarAngle - The desired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
% Outputs:
%   rgcDensityFit - handle of a fitting function that returns rgc density
%       values for the specified meridian angle across eccentricity in
%       polarAngle
%
% Options:
%   meridiansForKnotDefinition - We find that the knots must be defined
%       using only the 180 and 270 degree meridians to produce knots that
%       work for all meridians. This is a hack.

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@iscell);
p.addParameter('meridiansForKnotDefinition',[3,4],@isnumeric);
p.addParameter('splineKnots',20,@isnumeric);
p.addParameter('splineOrder',4,@isnumeric);

% Optional display params
p.addParameter('verbose',true,@islogical);
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(polarAngle,varargin{:})

%% sanity check the input
if p.Results.meridianAngles ~= 0
    error('This routine assumes that the first meridian is polar angle = 0');
end

% Loop across the cardinal meridians and combine the RGC data. Use this
% aggregate data to define the knots for the spline fit.
aggregatePosition=[];
aggregateDensity=[];
for mm=p.Results.meridiansForKnotDefinition(1):p.Results.meridiansForKnotDefinition(2)
    % load the empirical RGC density measured by Curcio
    [rgcDensitySqDeg, rgcNativeSupportPosDeg] = getCurcioRGCDensityByEccen(p.Results.meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(rgcDensitySqDeg));
    rgcNativeSupportPosDeg = rgcNativeSupportPosDeg(isvalididx);
    rgcDensitySqDeg = rgcDensitySqDeg(isvalididx);
    aggregatePosition=[aggregatePosition rgcNativeSupportPosDeg'];
    aggregateDensity = [aggregateDensity rgcDensitySqDeg'];
end

% perform a least-squares spline fit with the specified knots and order
BformSplineFit=spap2(p.Results.splineKnots, p.Results.splineOrder, aggregatePosition', aggregateDensity');
knots = BformSplineFit.knots;

% Loop across the cardinal meridians again, and now perform the spline fit
% with the specified knots
for mm=1:length(p.Results.meridianAngles)
    % load the empirical rgc density measured by Curcio
    [rgcDensitySqDeg, rgcNativeSupportPosDeg] = getCurcioRGCDensityByEccen(p.Results.meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(rgcDensitySqDeg));
    rgcNativeSupportPosDeg = rgcNativeSupportPosDeg(isvalididx);
    rgcDensitySqDeg = rgcDensitySqDeg(isvalididx);
    % Perform the spline fit
    BformSplineFit=spap2(knots, p.Results.splineOrder, rgcNativeSupportPosDeg', rgcDensitySqDeg');
    % convert from B-form to piecewise polynomial form
    ppFormSplineFits{mm} = fn2fm(BformSplineFit,'pp');
    % While we are here, record what is the smallest eccentricity for which
    % there is a non-zero RGC density measurement
    minEccNonZeroRGCDensity(mm) = rgcNativeSupportPosDeg(find(rgcDensitySqDeg ~= 0,1));
end

% We will now assemble a new coefficient matrix for the spline fit that is
% an interpolation between the coefficient matrices obtained for the
% cardinal meridians. This is done by going through each element of the
% coefficient matrix and obtaining the four values, one from each meridian.
% The value from polar angle 0 is replicated and assigned a polar angle of
% 360; this allows the interpolation to wrap-around. We then fit a 'pchip'
% (Piecewise Cubic Hermite Interpolating Polynomial) to the set of 5 values
% and interpolate to the called-for, intermediate polar angle.

% make an empty coefficient matrix
interpCoefs = zeros(size(ppFormSplineFits{1}.coefs));

% replicate the values for polar angle zero at polar angle 360 to allow a
% wrap-around fit
angleBase=[p.Results.meridianAngles 360];

% loop over each element of the coefficient matrix
for cc=1:numel(interpCoefs)
    % obtain the set of coefficients from across the meridians
    thisCoefByMeridian = cellfun(@(x) x.coefs(cc), ppFormSplineFits);
    % duplicate the first coefficient (assuming it is polar angle zero)
    thisCoefByMeridian = [thisCoefByMeridian thisCoefByMeridian(1)];
    % fit the pchip model
    coefByAngleFit = fit(angleBase',thisCoefByMeridian','pchipinterp');
    % evaluate the fitted model at the called for polar angle and place the
    % returned, interpolated coefficient value into the matrix we are
    % building
    interpCoefs(cc) = coefByAngleFit(polarAngle);
end

% Assemble our interpolated piecewise polynomial fit structure
ppFormSplineInterp = ppFormSplineFits{1};
ppFormSplineInterp.coefs = interpCoefs;

% Identify the smallest eccentricity value for which there was a non-zero
% RGC density in the empirical data in the closest cardinal meridian
[~,closestMeridian]=min(abs(p.Results.meridianAngles-polarAngle));


% Create an anonymous function using the interpolated spline. This function
% will return rgc density as a function of eccentricity in degrees.
% The second portion of the function forces the RGC density value to be
% zero if the eccentricity is in the range of the fovea with zero RGCs
rgcDensityFit = @(supportPosDeg) ...
    fnval(ppFormSplineInterp,supportPosDeg') .* ...
    (supportPosDeg'>=minEccNonZeroRGCDensity(closestMeridian));

end % function

