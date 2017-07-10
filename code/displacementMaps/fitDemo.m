%demo of fit functions

%% Fit Functions
%input angle
angle = 0;
% fit the RGC density
[ecc_mm,outParams_RGC,RGCdensityFit, scaleData] = fitRGCdensityDev(angle);
% fit the RF density -- need to convert mm to deg
RFfit = fitRFdensity(convert_mm_to_deg(ecc_mm),angle,scaleData);

%% Plot 
figure;
plot(ecc_mm,RGCdensityFit)
hold on
plot(ecc_mm,RFfit(convert_mm_to_deg(ecc_mm)))

%% Find K offset for RGC fit integral 
% create function 
shape    = outParams_RGC(1);
scale    = outParams_RGC(2);
location = outParams_RGC(3);

RGF_function = @(x1)(exp(-((x1-location)/scale).^-shape));
K_RGC = 0 - RGF_function(0);  

figure
plot(ecc_mm,(RGF_function(ecc_mm)+K_RGC)*scaleData)

%% Find K offset for RF fit integral 
a = RFfit.a;
b = RFfit.b;
c = RFfit.c;
d = RFfit.d;

RF_Function = @(x2)((a.*exp(b.*x2)./b) + (c.*exp(d.*x2)./d));
K_RF = (2*(14804.6)/scaleData) - RF_Function(0);

figure
plot(ecc_mm,(RF_Function(convert_mm_to_deg(ecc_mm))+K_RF)*scaleData);

%% Calculate Displacement 


