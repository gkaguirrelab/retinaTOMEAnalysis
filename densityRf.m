function Df = densityRf(r,theta)
% Density of recetive feilds
Rv0 = 0.011785;
Ro0 = 0.008333;

E20 = 20;

if theta >= 0 & theta <= 90
    E2v = interp1([0,90],[2.19,1.98],theta);
    
elseif theta > 90 & theta <= 180
    E2v = interp1([90,180],[1.98,2.26],theta);
   
elseif theta > 180 & theta <= 270
    E2v = interp1([180,270],[2.26,1.5],theta);
    
elseif theta > 270 & theta <=360
    E2v = interp1([270,360],[1.5,2.19],theta);
    
end

Rve = Rv0.*(1+r./E2v);
Roe =Ro0.*(1+r./E20);

k=1+(1.004 - 0.007209.*r + 0.001694.*r.^2 - 0.00003765.*r.^3).^-2;

Df = ((1.12+0.0273.*r).*k)./(1.155.*(Rve.^2-Roe.^2));

end