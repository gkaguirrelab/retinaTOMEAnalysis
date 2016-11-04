function Df = densityRf(r,theta,radDeg,smpPerDeg)
% Density of recetive feilds
%
% TO DO: add number of samples (per deg? or per 

Rv0 = 0.011785;
Ro0 = 0.008333;

E20 = 20;
a = 1;
for i = 1:size(theta,1);
    for ii = 1:size(theta,2)
        if theta(i,ii) >= 0 & theta(i,ii) <= 90
            E2v(i,ii) = interp1([0,90],[2.19,1.98],theta(i,ii));
           
        elseif theta(i,ii) > 90 & theta(i,ii) <= 180
            E2v(i,ii) = interp1([90,180],[1.98,2.26],theta(i,ii));
           
        elseif theta(i,ii) > 180 & theta(i,ii) <= 270
            E2v(i,ii) = interp1([180,270],[2.26,1.5],theta(i,ii));
             
        elseif theta(i,ii) > 270 & theta(i,ii) <=360
            E2v(i,ii) = interp1([270,360],[1.5,2.19],theta(i,ii));
        end
    end
end
Rve = Rv0.*(1+r./E2v);
Roe =Ro0.*(1+r./E20);

k=1+(1.004 - 0.007209.*r + 0.001694.*r.^2 - 0.00003765.*r.^3).^-2;

Df = ((1.12+0.0273.*r).*k)./(1.155.*(Rve.^2-Roe.^2));

lin = linspace(-radDeg,radDeg,2.*radDeg*smpPerDeg);
[x,y] = meshgrid(lin,lin);
mask=((x.^2+y.^2).^.5<=radDeg);

Df(mask~=1) = nan;

end