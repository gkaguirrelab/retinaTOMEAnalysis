%demo of fit functions

%% Fit Functions
%input angle
angle = 0;
% fit the RGC density
[ecc_mm,RGCfit] = fitRGCdensity(angle,6);
% fit the RF density -- need to convert mm to deg
RFfit = fitRFdensity(convert_mm_to_deg(ecc_mm),angle);

%% Plot 
figure;
plot(ecc_mm,exp(RGCfit(log(ecc_mm))))
hold on
plot(ecc_mm,RFfit(convert_mm_to_deg(ecc_mm)))

%% Sum and plot

RGCplot = exp(RGCfit(log(ecc_mm)));
RGCplot(1) = 0;
sum_RGCfit = cumsum(RGCplot);

sum_RFfit = cumsum(RFfit(convert_mm_to_deg(ecc_mm)));

figure;
plot(convert_mm_to_deg(ecc_mm),sum_RGCfit)
hold on
plot(convert_mm_to_deg(ecc_mm),sum_RFfit)

