function [fitParams, fitDisplacement] = fitExpExpToDisplacement(radii_mm, displacement, weights, varargin)
% function [fitParams, fitDisplacement] = fitExpExpToDisplacement(radii_mm, displacement, weights, varargin)
%
%  Fits an exponentiated exponential to the retinal ganglion cell displacement function.
%
%  http://home.iitk.ac.in/~kundu/paper101.pdf
%
%  Inputs:
%    radii_mm: vector with the eccentricity base of the data
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
p.addParameter('initialParams',[100,2],@isnumeric);
p.addParameter('displayPlot','none',@ischar);
p.parse(radii_mm, displacement, weights, varargin{:});

% Unpack the arguments
initialParams=p.Results.initialParams;
displayPlot=p.Results.displayPlot;

% Check that all passed vectors are either row or column vectors and the
% same length
if ~all(size(radii_mm)==size(displacement)) || ~all(size(radii_mm)==size(weights))
    error('The passed vectors must be the same length and same row/column order');
end

% Run the search
fitParams = fminsearch( (@(p) expExpModelFit(p, radii_mm, displacement, weights)), initialParams);
% Obtain the Gamma model fit at the passed eccentricities.
fitDisplacement = (fitParams(1).*fitParams(2)).*((1-exp(-1.*fitParams(1).*radii_mm)).^fitParams(2)).*exp(-1.*fitParams(1).*radii_mm);

% If the user requested a plot, give it to them
if strcmp(displayPlot,'full')
    figure;
    % Plot the data
    r1 = plot(radii_mm, displacement, '.r'); hold on;
    r2 = plot(radii_mm, fitDisplacement);
    % Make the plot pretty
    xlabel('Eccentricity (mm)');
    ylabel('Displaement (mm)');
end

end % main function


function E = expExpModelFit(params, radii_mm, displacement, weights)

% Error function, calculating the sum-of-squares for the data vs. the fit.
yhat = (params(1).*params(2)).*((1-exp(-1.*params(1).*radii_mm)).^params(2)).*exp(-1.*params(1).*radii_mm);

% Calculate the sums-of-square error
errorPreSum =((displacement - yhat).^2);
errorPreSum = errorPreSum .* weights;
E = sum(errorPreSum);

% Determine if the derivative of the gamma function contains points with a
% slope greater than unity. If so, inflate the error term by 10.
% This is to avoid a non-physiologic displacement of an RGC back
% past a more central, neighboring RGC.

firstDeriv=diff(yhat)./diff(radii_mm);
if max(firstDeriv) >=1
    E=E*1;
end


end % error calculation for gamma pdf model
