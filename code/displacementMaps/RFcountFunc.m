


nasal          = 2*(14804.6) .* (0.9729*((1+ecc_deg./1.084)).^-2)+(1-0.9729).*exp(-1.*ecc_deg./7.633);
%nasal          = (nasal .* (0.8928*((1+ecc_deg./41.03).^-1)))./scaleData;
nasal          = (nasal .* (0.8928*((1+ecc_deg./41.03).^-1)));
superior       = 2*(14804.6) * ( 0.9935*(1+ecc_deg/(1.035)).^-2+(1-0.9935)*exp(-1*ecc_deg/16.35));
%superior       = (superior .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;
superior       = (superior .* (0.8928*(1+ecc_deg./41.03).^-1));
temporal       = 2*(14804.6) * ( 0.9851*(1+ecc_deg/(1.058)).^-2+(1-0.9851)*exp(-1*ecc_deg/22.14));
%temporal       = (temporal .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;
temporal       = (temporal .* (0.8928*(1+ecc_deg./41.03).^-1));

inferior       = 2*(14804.6) * ( 0.996*(1+ecc_deg/(0.9932)).^-2+(1-0.996)*exp(-1*ecc_deg/12.13));
%inferior       = (inferior .* (0.8928*(1+ecc_deg./41.03).^-1))./scaleData;
inferior       = (inferior .* (0.8928*(1+ecc_deg./41.03).^-1));

curve_nasal    = fit(ecc_deg,nasal,'exp2','Exclude', find(isnan(nasal)));
curve_superior = fit(ecc_deg,superior,'exp2','Exclude', find(isnan(superior)));
curve_temporal = fit(ecc_deg,temporal,'exp2','Exclude', find(isnan(temporal)));
curve_inferior = fit(ecc_deg,inferior,'exp2','Exclude', find(isnan(inferior)));

% Find K offset for RF fit integral 
a = curve_superior.a;
b = curve_superior.b.*-1;
c = curve_superior.c;
d = curve_superior.d.*-1;

y = @(x)(2.*pi.*((x+1).*(-a.*exp(-b.*(x+1))/b - c.*exp(-d.*(x+1))/d) - (exp(-b.*(x+1))/(b.^2) + c.*exp(-d.*(x+1))/(d.^2))) - ...
    2.*pi.*(x.*(-a.*exp(-b.*x)/b - c.*exp(-d.*x)/d) - (exp(-b.*x)/(b.^2) + c.*exp(-d.*x)/(d.^2))));


%
figure
subplot(1,2,1)
plot(1:0.2:20,y(1:0.2:20).*scaleData)
title('RF Count')
xlabel('eccentricity (deg)')
ylabel('RF Count (receptive fields)')
subplot(1,2,2)
plot(1:0.2:20,cumsum(y(1:0.2:20).*scaleData))
title('Cumulative RF Count')
xlabel('eccentricity (deg)')
ylabel('Cumulative RF Count (receptive fields)')