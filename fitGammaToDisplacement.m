function [fitParams] = fitGammaToDisplacement(radii_mm, Displacement, weights)


%



%% define default parameters

% Initial guess for the parameters to define the Watson fit
%   lag
%   amplitude -- overall amplitude
%   gamma1 -- positive gamma parameter (roughly, time-to-peak in seconds)
%   gamma2 -- negative gamma parameter (roughly, time-to-peak in seconds)
%   Scale  -- scaling factor between the positive and negative gamma componenets

% parameters of the double-gamma hemodynamic filter (HRF)

initialParams(1) = 1;
initialParams(2) = 3;

% Run the optimization. We will return these params.
fitParams = fminsearch( (@(p) gammaModelFit(p, radii_mm, Displacement, weights)), initialParams);

% Obtain the Gamma model fit at the passed eccentricities. This might
% be used for plotting.

gammaFitToData = gammaModel(radii_mm, fitParams);

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
yhat = gammaModel(t, params);

% Calculate the sums-of-square error
errorPreSum =((y - yhat).^2);
errorPreSum = errorPreSum * weights;
E = sum(weights);

function H = gammaModel(t, params)

% Un-pack the passed parameters
params_amplitude = params(1);
params_gammaShape = params(2); 

% Generate the model. We return H.
H =  gampdf(t, params_gammaShape, params_amplitude);

