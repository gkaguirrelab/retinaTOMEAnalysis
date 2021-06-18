function [p, Yfit, fVal] = fitDensitySurface(Y,w,preFitAvgEccen,simplePolarModel,p0,supportDeg,maxSupportDeg,refEccen,refDensity)
% Fit a multi-parameter surface to cone density data
%
% Syntax:
%   [p, Yfit, fVal] = fitDensitySurface(Y,w,preFitAvgEccen,simplePolarModel,p0,supportDeg,maxSupportDeg)
%
% Description:
%   A model of cone density across eccentricity and polar angle.
%   Variationin density across eccentricity is modeled as the sum of two
%   exponentials, and thus four parameters. Variation across polar angle is
%   modeled as a multiplicative adjustment of density as a function of
%   polar angle under the control of Fourier functions. These
%   functions include:
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
%

arguments
    Y (:,:) {mustBeNumeric}
    w (:,:) {mustBeNumeric} = ones(size(Y))
    preFitAvgEccen (1,1) = true;
    simplePolarModel (1,1) = true;
    p0 (1,20) {mustBeNumeric} = [1.4484e+03, -0.0825, 9.3044e+03, -1.3022, 31.8197, 0.0617, 8.9892, 0.1614, 0, 0.0601, 5.2897, 1.0714, -9.2917, 0.0958, 4.0656, 0.4412, -3.3665, 0.1081, 2.8328, 1.999]
    supportDeg (1,:) {mustBeNumeric} = 0:0.0078:0.0078*(size(Y,1)-1)
    maxSupportDeg (1,1) {mustBeNumeric} = 15
    refEccen (1,1) = 30 
    refDensity (1,1) = 500
end

%% pBlock and mBlock settings
% These are the parameters of the model. The pBlock defines the form of the
% sum of exponential fit. The mBlock defines the parameters of a gamma
% function modulation of a sinusoidal variation in density across polar
% angle.

pBlockLB = [0,-5,0,-5];
pBlockUB = [5e4,0,5e4,0];

mBlockLB = [-35 0 2 0.01];
mBlockUB = [35 1 20 2];

% The number of Fourier components that models variation across polar angle
nFourier = 4;

% The length of the support, pulled out here for code clarity
supportLength = length(supportDeg);


%% Initial fit to mean density profile

if preFitAvgEccen

    % Density, eccentricity, polar angle
    Ym = nanmean(Y,1);
    wm = nanmean(w,1);
    X = supportDeg;
    P = zeros(size(X))+90;
    
    % Bounds
    lb = [pBlockLB, p0(5:end)];
    ub = [pBlockUB, p0(5:end)];

    % Options
    options = optimoptions('fmincon','Display','off','Algorithm','interior-point');

    % objective
    validIdx = ~isnan(Ym);
    myObj = @(p) norm( wm(validIdx).* (Ym(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );
    
    % search
    pm = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % Update p0
    p0(1:4) = pm(1:4);

end


%% Fit polar angle variation

% Density, eccentricity, polar angle
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);
validIdx = ~isnan(Y);

if simplePolarModel
    
    % p0 and bounds
    lb = [pBlockLB, -30 0];
    ub = [pBlockUB, 30 2];
    
    % assemble p function
    pFull = @(p) [p(1:4) p0(5)+p(5) p0(6)*p(6) p0(7:8) p0(9)+p(5) p0(10)*p(6) p0(11:12) p0(13)+p(5) p0(14)*p(6) p0(15:16) p0(17)+p(5) p0(18)*p(6) p0(19:20)];
    
    % truncate p0
    p0 = [p0(1:4) 0 1];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportDeg,pFull(p))) );
    
else
    
    % p0 and bounds
    lb = [pBlockLB, repmat(mBlockLB,1,nFourier)];
    ub = [pBlockUB, repmat(mBlockUB,1,nFourier)];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - coneDensityModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );
    
end

% non-linear constraint for the asymptotic cone density
myNonlcon = @(p) asymptoteDensity(p,refEccen,refDensity);

% search
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','UseParallel',true);
[p, fVal] = fmincon(myObj,p0,[],[],[],[],lb,ub,myNonlcon,options);

if simplePolarModel
    p = pFull(p);
end

% generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:)=coneDensityModel(X,P,maxSupportDeg,p);

end


function [c,ceq] = asymptoteDensity(p,refEccen,refDensity)

% Decompose p into individual variables
a = p(1);   % scale first exponential
b = p(2);   % time constant first exponential
c = p(3);   % scale second exponential
d = p(4);   % time constant second exponential

% The modeled density at the reference eccentricity
density = (a.*exp(b.*refEccen)+c.*exp(d.*refEccen));

% Constraint to keep density close to the reference density
ceq = density - refDensity;
c = [];

end