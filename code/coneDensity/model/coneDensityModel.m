function density = coneDensityModel(x,angle,maxX,p)
% A parameterized model of retinal cone density as a function of eccen / PA
%
% Syntax:
%   density = coneDensityModel(x,angle,maxX,p)
%
% Description:
%   A model of cone density across eccentricity and polar angle.
%   Variation in density across eccentricity is modeled as the sum of two
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
%   x                     - Vector. Retinal eccentricity (mm or degrees)
%   angle                 - Vector. Retinal polar angle in degrees
%   maxX                  - Scalar. The reference eccentricity location for
%                           setting up the gamma functions. This may be a
%                           value that is larger or smaller than the
%                           maximum value in supportDeg.
%   p                     - 1x20 vector of params of the model.
%
% Outputs:
%   density               - Vector. Cone density.
%

% Decompose p into individual variables
a = p(1);   % scale first exponential
b = p(2);   % time constant first exponential
c = p(3);   % scale second exponential
d = p(4);   % time constant second exponential
ph1 = p(5); % angle adjustment of cos1
f11 = p(6); % gamma amplitude for cos1
f12 = p(7); % gamma shape for cos1
f13 = p(8); % gamma scale for cos1
ph2 = p(9); % angle adjustment of sin1 
f21 = p(10);% gamma amplitude for sin1
f22 = p(11);% gamma shape for sin1
f23 = p(12);% gamma scale for sin1
ph3 = p(13);% angle adjustment of cos2
f31 = p(14);% gamma amplitude for cos2
f32 = p(15);% gamma shape for cos2
f33 = p(16);% gamma scale for cos2
ph4 = p(17);% angle adjustment of cos4
f41 = p(18);% gamma amplitude for cos4
f42 = p(19);% gamma shape for cos4
f43 = p(20);% gamma scale for cos4

% Assemble the gamma pdf functions, one for each Fourier modulation
g1 = f11.*gampdf(x,f12,f13)./max(gampdf(0:0.01:maxX,f12,f13));
g2 = f21.*gampdf(x,f22,f23)./max(gampdf(0:0.01:maxX,f22,f23));
g3 = f31.*gampdf(x,f32,f33)./max(gampdf(0:0.01:maxX,f32,f33));
g4 = f41.*gampdf(x,f42,f43)./max(gampdf(0:0.01:maxX,f42,f43));

% Assemble the modulatory component of the model. This is the
% multiplicative variation in cone density as a function of polar angle.
% The angle of the cos1 and sin1 components are shifted 180 degrees, so
% that the scaling parameter of the gamma functions will be in the positive
% domain for empirical data.
m = g1.*cosd(angle+ph1+180) + ...
    g2.*sind(angle+ph2+180) + ...
    g3.*cosd(angle.*2+ph3) + ...
    g4.*cosd(angle.*4+ph4);

% The final density model is comprised of a double exponential as a
% function of eccentricity, modulated by the polar angle.
density = (m+1).*(a.*exp(b.*x)+c.*exp(d.*x));

end

