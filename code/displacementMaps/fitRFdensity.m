function expFitOut = fitRFdensity(ecc_deg,angle,scaleData)
%fitRFdensity -- Estimates a RGC receptive field density function at a given angle on the retina.
%
% Description:
%   This function produces a function that esimates the retinal ganglion cell
%   receptive field density as a function of eccentricity in degrees .This
%   uses the Watson (2014) fits. This returns a function that estimates
%   this density at the desired input angle by taking a weighted average of
%   the fit parameters.
%
% Inputs:
%   ecc_deg = Sample positions along a meridian in degrees.
%   angle   = The dedsired angle of the density function on the retinal field.
%             (0=nasal;90=superior;180=temporal;270=inferior)
%
% Outputs:
%   expFitOut = Function that estimates the RGC recetive field denstity at
%               the input anlge as a funciton of eccentricity (deg).
%
% MAB 2017


% Use the equation in Wastson 2014 to estimate the ganglion cell receptive
% field density as a function of eccentricity in degrees.
% This estimate is then fit with a two term expontial of the form f(x) = a*exp(b*x) + c*exp(d*x)

% Generating mRGCf density from Watson 2014 scaled by max values from fitRGCdensityDev
nasal          = 2*(14804.6) .* (0.9729*((1+ecc_deg./1.084)).^-2)+(1-0.9729).*exp(-1.*ecc_deg./7.633);
nasal          = (nasal .* (0.8928*((1+ecc_deg./41.03).^-1)))./scaleData;

superior       = 2*(14804.6) * ( 0.9935*(1+ecc_deg/(1.035)).^-2+(1-0.9935)*exp(-1*ecc_deg/16.35));
superior       = (superior .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;

temporal       = 2*(14804.6) * ( 0.9851*(1+ecc_deg/(1.058)).^-2+(1-0.9851)*exp(-1*ecc_deg/22.14));
temporal       = (temporal .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;

inferior       = 2*(14804.6) * ( 0.996*(1+ecc_deg/(0.9932)).^-2+(1-0.996)*exp(-1*ecc_deg/12.13));
inferior       = (inferior .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;

% Fit with a two term expontial of the form f(x) = a*exp(b*x) + c*exp(d*x)
% to get the parameters a,b,c,d
curve_nasal    = fit(ecc_deg,nasal,'exp2','Exclude', find(isnan(nasal)));
curve_superior = fit(ecc_deg,superior,'exp2','Exclude', find(isnan(superior)));
curve_temporal = fit(ecc_deg,temporal,'exp2','Exclude', find(isnan(temporal)));
curve_inferior = fit(ecc_deg,inferior,'exp2','Exclude', find(isnan(inferior)));

% Set an output function handle with the purpose of overwriting the parameters.
expFitOut = curve_nasal;

% Take a weighed average of the parameters of the fit exp2. The
% weights are thet fraction of the input angle for both meridians that flank the
% angle of interest.
if angle >= 0 && angle < 90;
    superiorFrac = angle/90;
    nasalFrac  = 1 - superiorFrac;
    for i = ['a','b','c','d']
        eval(sprintf('expFitOut.%s = nasalFrac.*curve_nasal.%s + superiorFrac.*curve_superior.%s;',i,i,i))
    end
elseif angle >= 90 && angle < 180;
    temporalFrac = (angle-90)/90;
    superiorFrac  = 1 - temporalFrac;
    for i = ['a','b','c','d']
        eval(sprintf('expFitOut.%s = superiorFrac.*curve_superior.%s + temporalFrac.*curve_temporal.%s;',i,i,i))
    end
elseif angle >= 180 && angle < 270;
    inferiorFrac = (angle-180)/90;
    temporalFrac= 1 - inferiorFrac;
    for i = ['a','b','c','d']
        eval(sprintf('expFitOut.%s = temporalFrac.*curve_temporal.%s + inferiorFrac.*curve_inferior.%s;',i,i,i))
    end
elseif angle >= 270 && angle < 360;
    nasalFrac = (angle-270)/90;
    inferiorFrac= 1 - nasalFrac;
    for i = ['a','b','c','d']
        eval(sprintf('expFitOut.%s = inferiorFrac.*curve_inferior.%s + nasalFrac.*curve_nasal.%s;',i,i,i))
        
    end
end


end