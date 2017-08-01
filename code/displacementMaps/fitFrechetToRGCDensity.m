function [fitParams, RGCdensityFit] = fitFrechetToRGCDensity(ecc_mm, rgcDensity, weights, varargin)
% function [fitParams, fitDisplacement] = fitGammaToDisplacement(radii_mm, displacement, weights, varargin)
%
%  Fits a gamma pdf to the retinal ganglion cell displacement function. The
%  use of a gamma pdf for the fit is driven both by the appearance of the
%  data and the prior use of this function in Watson 2014 JoV, Eq 5.
%
%  Inputs:
%    ecc_mm: vector with the eccentricity base of the data
%    displacement: vector with the RGC displacement from the fovea
%    weights: vector that provides a weighting function for the error
%      term as a function of eccentricity. The use of a weighting function
%      is motivated by the observation that the Turpin...McKendrick
%      solution that we implement for the calculation of displacement does
%      not always yield a measured displacement function that reaches zero.
%      For this reason we fit a gamma function to the upslope of the
%      displacement function, and discount the more eccentric values.
%
%  Optional inputs:
%    initialParams: a key value pair with a two element vector (shape,
%       scale)
%    displayPlot: a key-value pair with values of:
%       full - produce plots of the fit
%       none - (default) no plots
%
%  Outputs
%    fitParams: the parameters of the best fit gamma pdf (shape, scale)
%    fitDisplacement: the gamma function itself
%
%
%  Demo - create a gamma pdf (plus noise), a uniform weighting
%    function, and then fit it
%
%    radii_mm=0:0.005:5;
%    displacement=gampdf(radii_mm,3,.2)+(rand(1,1001)-0.5)/5;
%    weights=radii_mm.*0+1;
%
%    [fitParams, fitDisplacement] = fitGammaToDisplacement(radii_mm, ...
%      displacement, weights, 'displayPlot', 'full');

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('radii_mm',@isnumeric);
p.addRequired('displacement',@isnumeric);
p.addRequired('weights',@isnumeric);
p.addParameter('initialParams',[1,1,0],@isnumeric);
p.addParameter('displayPlot','none',@ischar);
p.parse(ecc_mm, rgcDensity, weights, varargin{:});

% Unpack the arguments
initialParams=p.Results.initialParams;
displayPlot=p.Results.displayPlot;

% Check that all passed vectors are either row or column vectors and the
% same length
if ~all(size(ecc_mm)==size(rgcDensity)) || ~all(size(ecc_mm)==size(weights))
    error('The passed vectors must be the same length and same row/column order');
end

options=optimset('fminsearch');
options=optimset(options,'Display','none');

% Run the search
fitParams = fminsearch( (@(p) modelFitError(p, ecc_mm, rgcDensity, weights)), initialParams, options);

% Obtain the Frechnet model fit at the passed eccentricities.
RGCdensityFit = frechnetPDF(ecc_mm,fitParams(1),fitParams(2),fitParams(3));

% If the user requested a plot, give it to them
if strcmp(displayPlot,'full')
    figure;
    % Plot the data
    r1 = plot(ecc_mm, rgcDensity, '.r'); hold on;
    r2 = plot(ecc_mm, RGCdensityFit);
    % Make the plot pretty
    xlabel('Eccentricity (mm)');
    ylabel('Density (counts/mm2)');
end

end % main function


function y = frechnetPDF( x, shape, scale, location )

y = (shape./scale)*((x-location)./scale).^(-1-shape).* exp( -((x-location)./scale).^(-shape));


end


function E = modelFitError(params, radii_mm, rgcDensity, weights)

% Error function, calculating the sum-of-squares for the data vs. the fit.
yhat = frechnetPDF(radii_mm,params(1),params(2),params(3));

% Calculate the sums-of-square error
errorPreSum =((rgcDensity - yhat).^2);
errorPreSum = errorPreSum .* weights;
E = nansum(errorPreSum);

end
