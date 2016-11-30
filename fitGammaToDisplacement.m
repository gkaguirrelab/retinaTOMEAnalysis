function [fitParams] = fitGammaToDisplacement(radii_mm, Displacement, weights)

% parameters of a gamma function (shape and scale)

initialParams(1) = 3; % shape
initialParams(2) = 1; % scale

% Run the optimization. We will return these params.
fitParams = fminsearch( (@(p) gammaModelFit(p, radii_mm, Displacement, weights)), initialParams);

% Obtain the Gamma model fit at the passed eccentricities. This might
% be used for plotting.

gammaFitToData = gampdf(radii_mm,fitParams(1),fitParams(2));

% If the user requested a plot, give it to them

displayFitPlot=true;
if displayFitPlot
    figure;
    % Plot the data
    r1 = plot(radii_mm, Displacement, 'sr', 'MarkerFaceColor', 'r'); hold on;
    r2 = plot(radii_mm, gammaFitToData);
    % Make the plot pretty
    xlabel('Eccentricity (mm)');
    ylabel('Displaement (mm)');
end




function E = gammaModelFit(params, t, y, weights)

% Error function, calculating the sum-of-squares for the data vs. the fit.
yhat = gampdf(t,params(1),params(2));

% Calculate the sums-of-square error
errorPreSum =((y - yhat).^2);
errorPreSum = errorPreSum .* weights;
E = sum(errorPreSum);
