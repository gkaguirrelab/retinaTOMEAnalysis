%demo of fit functions
clear all
%% Fit Functions
%input angle 
angle = 0;
% fit the RGC density
[ecc_deg,outParams_RGC,RGCdensityFit, scaleData] = fitRGCdensityDev(angle);
% fit the RF density -- need to convert mm to deg
RFfit = fitRFdensity(ecc_deg,angle,scaleData);

%% Plot 
figure;
plot(ecc_deg,RGCdensityFit)
hold on
plot(ecc_deg,RFfit(ecc_deg))

%% Find K offset for RGC fit integral 
% create function 
shape    = outParams_RGC(1);
scale    = outParams_RGC(2);
location = outParams_RGC(3);
RGC_function = @(x1)(exp(-((x1-location)/scale).^-shape));
K_RGC = 0 - RGC_function(0);  




%% Find K offset for RF fit integral 
a = RFfit.a;
b = RFfit.b;
c = RFfit.c;
d = RFfit.d;

RF_Function = @(x2)((a.*exp(b.*x2)./b) + (c.*exp(d.*x2)./d));
K_RF = (2*(14804.6)./scaleData) - RF_Function(0);


%% plot
figure
plot(ecc_deg,(RF_Function(convert_mm_to_deg(ecc_deg))+K_RF),'r');
hold on 
plot(ecc_deg,((RGC_function(ecc_deg)+K_RGC)),'b')


%% Calculate Displacement 

% 
Drf = RF_Function(convert_mm_to_deg(ecc_deg))+K_RF;

RGC_postition = (-log(Drf-K_RGC)).^(-1./shape).*(location.*(-log(Drf-K_RGC)).^(1./shape) + scale);
 
displacement =RGC_postition - ecc_deg;

figure 
plot(ecc_deg,displacement)

