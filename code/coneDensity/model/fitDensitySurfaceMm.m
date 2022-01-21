function [p, Yfit, fVal, RSquared, nonlcon, polarTheta, polarMultiplier] = fitDensitySurfaceMm(Y,w,preFitAvgEccen,simplePolarModel,usePeripheralDensityConstraint,useFovealDensityConstraint,p0,supportMm,maxSupportMm,refPeripheralLocation,refPeripheralDensity,refFovealDensity)
% Fit a multi-parameter surface to cone density data
%
% Syntax:
%   [p, Yfit, fVal, RSquared, nonlcon, polarTheta, polarMultiplier] = fitDensitySurfaceMm(Y,w,preFitAvgEccen,simplePolarModel,usePeripheralDensityConstraint,useFovealDensityConstraint,p0,supportMm,maxSupportMm,refPeripheralLocation,refPeripheralDensity,refFovealDensity)
%
% Description:
%   A model of cone density across eccentricity and polar angle.
%   Variation in density across eccentricity is modeled as the sum of two
%   exponentials, and thus four parameters. Variation across polar angle is
%   modeled as a multiplicative adjustment of density as a function of
%   polar angle under the control of Fourier functions. These functions
%   include:
%     - a sine and cosine at the fundamental frequency, which model
%       variation in density between the superior and inferior retina, and
%       betwween the nasal and temporal retina.
%     - a cosine at the second frequency, which models variation in density
%       between the vertical and horizontal meridians
%     - a cosine at the fourth frequency, which models variation in density
%       between the meridian and non-meridian angles.
%
%   The multiplicative effect of each of these Fourier components is
%   subject to a magnitude scaler, which itself varies as a function of
%   eccentricity under the control of a two-parameter gamma function.
%   Finally, each of the Fourier components may be independently phase
%   shifted, with the exception of first sine and cosine, which are
%   assigned the same phase shift.
%
% Inputs:
%   Y                     - An n x n polar matrix, with the rows being
%                           eccentricity and the columns polar angle.
%   w                     - An n x n matrix of weights for the fit.
%   preFitAvgEccen        - Logical. Determines if the fit is initialized
%                           by first fitting an exponential to the mean
%                           density as a function of eccentricity,
%                           collapsing across polar angle.
%   simplePolarModel      - Logical. When set to true, the parameters of
%                           the polar angle variation in density are held
%                           constant, subject only to an overall
%                           multiplicative scaling of the strength of the
%                           modulation, and an angle that can rotate the
%                           entire polar model.
%   p0                    - 1x20 vector. The initial values of the model.
%   supportDeg            - 1xn vector. The support for the data in the 
%                           eccentricity direction, in units of degrees.
%   maxSupportDeg         - Scalar. The reference eccentricity location for
%                           setting up the gamma functions. This may be a
%                           value that is larger or smaller than the
%                           maximum value in supportDeg.
%
% Outputs:
%   p                     - 1x20 vector. The parameters of the model fit.
%   Yfit                  - nxn matrix. The model fit
%   fVal                  - Scalar. The fit error.
%   RSquared              - 1x4 vector. A set of R-squared values for:
%                             - RSquaredFull
%                             - RSquaredExponentialOnly
%                             - RSquaredPolarOnly
%                             - RSquaredPolarResiduals
%   polarMultiplier       - Scalar. When simplePolarModel is set to true,
%                           this is the multiplicative scaler that is
%                           applied to the effect of polar angle upon the
%                           fit. If simplePolarModel is set to false, then
%                           this will have a value of nan.
%

arguments
    Y (:,:) {mustBeNumeric}
    w (:,:) {mustBeNumeric} = ones(size(Y))
    preFitAvgEccen (1,1) = true;
    simplePolarModel (1,1) = true;
    usePeripheralDensityConstraint (1,1) = false;
    useFovealDensityConstraint (1,1) = false;
    p0 (1,20) {mustBeNumeric} = [ ...
        1.6349e+03,  -0.2816, 8.8862e+03, -5, ...
        53.6134    0.1128    4.6850    0.1085, ...
       -11.2582    0.0560    1.6059    1.0000, ...
      -10.8416    0.0832    2.2149    0.2996, ...
       -5.4307    0.0661    3.8713    0.1751]
    supportMm (1,:) {mustBeNumeric} = 0:0.0025:0.0025*(size(Y,1)-1)
    maxSupportMm (1,1) {mustBeNumeric} = 5
    refPeripheralLocation (1,1) = 5
    refPeripheralDensity (1,1) = 400
    refFovealDensity (1,1) = 14426; % Median peak cone density from Table S1 of Reiniger et al., 2021 Current Biology
end

%% pBlock and mBlock settings
% These are the parameters of the model. The pBlock defines the form of the
% sum of exponential fit. The mBlock defines the parameters of a gamma
% function modulation of a sinusoidal variation in density across polar
% angle.

pBlockLB = [0,-10,0,-10];
pBlockUB = [5e4,0,5e4,0];

mBlockLB = [-90  0.0  1.0  0.05];
mBlockUB = [ 90  0.2 15.0  1.5];


% The number of Fourier components that models variation across polar angle
nFourier = 4;

% The length of the support, pulled out here for code clarity
supportLength = length(supportMm);


%% Initial fit to mean density profile

if preFitAvgEccen

    % Density, eccentricity, polar angle
    Ym = nanmean(Y,1);
    wm = nanmean(w,1);
    X = supportMm;
    P = zeros(size(X))+90;
    
    % Bounds
    lb = [pBlockLB, p0(5:end)];
    ub = [pBlockUB, p0(5:end)];

    % Options
    options = optimoptions('fmincon','Display','off','Algorithm','interior-point');

    % objective
    validIdx = ~isnan(Ym);
    myObj = @(p) norm( wm(validIdx).* (Ym(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportMm,p)) );
    
    % search
    pm = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % Update p0
    p0(1:4) = pm(1:4);

end


%% Fit polar angle variation

% Density, eccentricity, polar angle
X = repmat(supportMm,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);
validIdx = ~isnan(Y);

if simplePolarModel
    
    % p0 and bounds
    lb = [pBlockLB, -30 -2];
    ub = [pBlockUB, 30 2];
    
    % assemble p function
    pFull = @(p) [p(1:4) p0(5)+p(5) p0(6)*p(6) p0(7:8) p0(9)+p(5) p0(10)*p(6) p0(11:12) p0(13)+p(5) p0(14)*p(6) p0(15:16) p0(17)+p(5) p0(18)*p(6) p0(19:20)];
    
    % truncate p0
    p0 = [p0(1:4) 0 1];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportMm,pFull(p))) );
    
else
    
    % p0 and bounds
    lb = [pBlockLB, repmat(mBlockLB,1,nFourier)];
    ub = [pBlockUB, repmat(mBlockUB,1,nFourier)];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportMm,p)) );
    
end

% non-linear constraint for the asymptotic cone density
myNonlcon = [];
if usePeripheralDensityConstraint && useFovealDensityConstraint
    myNonlcon = @(p) constrainFovealAndPeripheralDensity(p,refPeripheralLocation,refPeripheralDensity,refFovealDensity);
end
if usePeripheralDensityConstraint && ~useFovealDensityConstraint
    myNonlcon = @(p) constrainPeripheralDensity(p,refPeripheralLocation,refPeripheralDensity);
end

% search
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','UseParallel',true);
[p, fVal] = fmincon(myObj,p0,[],[],[],[],lb,ub,myNonlcon,options);

% Obtain the non-linear constraint value
if usePeripheralDensityConstraint
    nonlcon = myNonlcon(p);
else
    nonlcon = nan;
end

% Expand the p vector if needed
if simplePolarModel
    polarTheta = p(5);
    polarMultiplier = p(6);
    p = pFull(p);
else
    polarTheta = nan;
    polarMultiplier = nan;
end

% Generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:) = coneDensityModel(X,P,maxSupportMm,p);

% Calculate the R-squared value for the exponential and polar variation
% components
pReduced = p;
pReduced([6 10 14 18])=0;
YfitExponentialOnly = nan(supportLength,supportLength);
YfitExponentialOnly(:,:) = coneDensityModel(X,P,maxSupportMm,pReduced);
YfitPolarOnly = Yfit-YfitExponentialOnly;

% Only evaluate the correlation where there are data values
goodIdx = ~isnan(Y(:));

% Calculate a few metrics of the fit, including the R-squared value for a
% few different model components
RSquaredFull = corr(Yfit(goodIdx),Y(goodIdx))^2;
RSquaredExponentialOnly = corr(YfitExponentialOnly(goodIdx),Y(goodIdx))^2;
RSquaredPolarOnly = corr(YfitPolarOnly(goodIdx),Y(goodIdx))^2;

YPolarResiduals = Y - YfitExponentialOnly;
RSquaredPolarResiduals = corr(YfitPolarOnly(goodIdx),YPolarResiduals(goodIdx))^2;


RSquared = [RSquaredFull; RSquaredExponentialOnly; RSquaredPolarOnly; RSquaredPolarResiduals];

end


function [c,ceq] = constrainFovealAndPeripheralDensity(p,refPeripheralLocation,refPeripheralDensity,refFovealDensity)

% Decompose p into individual variables
a = p(1);   % scale first exponential
b = p(2);   % time constant first exponential
c = p(3);   % scale second exponential
d = p(4);   % time constant second exponential

% The modeled density at the peripheral reference eccentricity, and the
% modeled density at the foveal center
modelPeripheralDensity = (a.*exp(b.*refPeripheralLocation)+c.*exp(d.*refPeripheralLocation));
modelFovealDensity = (a.*exp(b.*1e-4)+c.*exp(d.*1e-4));

% Constraint to keep peripheral density above the reference peripheral
% density, and to match the peak density to the reference peak density
c = refPeripheralDensity - modelPeripheralDensity;
ceq = refFovealDensity - modelFovealDensity;

end



function [c,ceq] = constrainPeripheralDensity(p,refPeripheralLocation,refPeripheralDensity)

% Decompose p into individual variables
a = p(1);   % scale first exponential
b = p(2);   % time constant first exponential
c = p(3);   % scale second exponential
d = p(4);   % time constant second exponential

% The modeled density at the reference eccentricity
modelPeripheralDensity = (a.*exp(b.*refPeripheralLocation)+c.*exp(d.*refPeripheralLocation));

% Constraint to keep peripheral density above the reference peripheral
% density,
c = refPeripheralDensity - modelPeripheralDensity;
ceq = [];

end