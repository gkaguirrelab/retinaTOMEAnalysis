function [coneDensityFit] = getConeDensityFit(polarAngle, varargin)
% getConeDensityFit(angle)
%
% This routine returns a fit to cone density data. Raw cone density data
% are taken from Curcio et al (1990). Fits to the cardinal meridians are
% obtained. Interpolation over the parameters of the fit are used to
% produce a function that returns cone density for an arbitrary meridian
% angle.
%
% Inputs:
%   polarAngle - The desired angle of the density function on the retinal field.
%                (0=nasal;90=superior;180=temporal;270=inferior)
% Outputs:
%   coneDensityFit - handle of a fitting function that returns cone density
%       values for the specified meridian angle across eccentricity in
%       polarAngle
%
% Options:
%   splineOnly - if one of the cardinal meridians has been requested, then
%       setting this parameter to true will cause the returned fit function
%       to be to a spline fit of the cone density data.

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('polarAngle',@isnumeric);

% Optional analysis params
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@iscell);
p.addParameter('splineOnly',false,@islogical);

% Optional display params
p.addParameter('verbose',true,@islogical);
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(polarAngle,varargin{:})

% Define a fitType that is a modified power function
customFunc = fittype('a.*(x+b).^c+d.*x+e','independent','x','dependent','y');

% Loop across the cardinal meridians and fit the raw cone data with this
% function
for mm=1:length(p.Results.meridianAngles)
    % load the empirical cone density measured by Curcio
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(p.Results.meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    
    startPoint = [max(coneDensitySqDeg) 1 1 1 1];
    customFit = fit(coneNativeSupportPosDeg', coneDensitySqDeg', customFunc, 'StartPoint',startPoint);
    paramLabels = coeffnames(customFit);
    meridianFitParams(mm,:)=cellfun(@(x) customFit.(x), paramLabels);
end

% Loop through the fit parameters. For each, interpolate to the requested
% polarAngle using a pchipinterp fit and interpolation

% replicate the values for polar angle zero at polar angle 360 to allow a
% wrap-around fit
meridianFitParams=[meridianFitParams;meridianFitParams(1,:)];
angleBase=[p.Results.meridianAngles 360];

for pp=1:length(paramLabels)
    paramFit = fit(angleBase',meridianFitParams(:,pp),'pchipinterp');
    coneDensityFitParams(pp)=paramFit(polarAngle);
end

if p.Results.splineOnly
    if sum([0 90 180 270]==polarAngle)~=1
        error('Spline fit only available for cardinal meridians');
    else
        % Create an fit function using a spline fit to the cone data
        coneDensityFit = fit(coneNativeSupportPosDeg',coneDensitySqDeg','smoothingspline','SmoothingParam', 1);
    end
else
    % Create an anonymous function using the interpolated params
    coneDensityFit = @(x) customFunc(coneDensityFitParams(1), coneDensityFitParams(2), coneDensityFitParams(3), coneDensityFitParams(4), coneDensityFitParams(5), x');
end

end % function

