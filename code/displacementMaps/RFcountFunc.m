function [CountPerRingRF] = RFcountFunc(RFfit,verbose)

% Find K offset for RF fit integral 
a = RFfit.a;
b = RFfit.b.*-1;
c = RFfit.c;
d = RFfit.d.*-1;

CountPerRingRF = @(x)(2.*pi.*((x+1).*(-a.*exp(-b.*(x+1))/b - c.*exp(-d.*(x+1))/d) - (exp(-b.*(x+1))/(b.^2) + c.*exp(-d.*(x+1))/(d.^2))) - ...
    2.*pi.*(x.*(-a.*exp(-b.*x)/b - c.*exp(-d.*x)/d) - (exp(-b.*x)/(b.^2) + c.*exp(-d.*x)/(d.^2))));


%
figure
subplot(1,2,1)
plot(1:0.2:20,CountPerRingRF(1:0.2:20).*scaleData)
title('RF Count')
xlabel('eccentricity (deg)')
ylabel('RF Count (receptive fields)')
subplot(1,2,2)
plot(1:0.2:20,cumsum(CountPerRingRF(1:0.2:20).*scaleData))
title('Cumulative RF Count')
xlabel('eccentricity (deg)')
ylabel('Cumulative RF Count (receptive fields)')