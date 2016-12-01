function [fitParams, fitDisplacement] = fitGammaToDisplacement(radii_mm, displacement, weights, varargin)
% function [fitParams, fitDisplacement] = fitGammaToDisplacement(radii_mm, displacement, weights, varargin)
%
%  Fits a gamma pdf to the retinal ganglion cell displacement function. The
%  use of a gamma pdf for the fit is driven both by the appearance of the
%  data and the prior use of this function in Watson 2014 JoV, Eq 5.
%
%  Inputs:
%    radii_mm: row vector with the eccentricity base of the data
%    displacement: row vector with the RGC displacement from the fovea
%    weights: row vector that provides a weighting function for the error
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

%% Parse vargin for options passed here
p = inputParser;
p.addRequired('radii_mm',@isnumeric);
p.addRequired('displacement',@isnumeric);
p.addRequired('weights',@isnumeric);
p.addParameter('initialParams',[3,1],@isnumeric);
p.addParameter('displayPlot','none',@ischar);
p.parse(thePacket,varargin{:});

% Unpack the arguments
initialParams=p.Result.initialParams;
displayPlot=p.Result.displayPlot;

% Run the optimization
fitParams = fminsearch( (@(p) gammaModelFit(p, radii_mm, displacement, weights)), initialParams);

% Obtain the Gamma model fit at the passed eccentricities.
fitDisplacement = gampdf(radii_mm,fitParams(1),fitParams(2));

% If the user requested a plot, give it to them
if strcmp(displayPlot,'full');
    figure;
    % Plot the data
    r1 = plot(radii_mm, displacement, 'sr', 'MarkerFaceColor', 'r'); hold on;
    r2 = plot(radii_mm, gammaFitToData);
    % Make the plot pretty
    xlabel('Eccentricity (mm)');
    ylabel('Displaement (mm)');
end

end % main function


function E = gammaModelFit(params, radii_mm, displacement, weights)

% Error function, calculating the sum-of-squares for the data vs. the fit.
yhat = gampdf(radii_mm,params(1),params(2));

% Calculate the sums-of-square error
errorPreSum =((displacement - yhat).^2);
errorPreSum = errorPreSum .* weights;
E = sum(errorPreSum);

end % error calculation for gamma pdf model
