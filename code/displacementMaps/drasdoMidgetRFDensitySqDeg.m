function RFdensitySqDeg = drasdoMidgetRFDensitySqDeg(supportPosDeg, angle)
% Genrates a vector of recetive feild dentity 
%
% From Eq. 7 in Drasdo et al. 2007 Vision Research.
% Formula is in Degrees 
% 
% The Drasdo equations are cast in terms of position in the visual field in
% degrees eccentricity and counts per square degree.
%
% The output of this routine will be in units of square degrees at an
% eccentricity given in degrees. The values, however, will be for that
% position in the retinal coordinate space, not the visual.
%
% The values returned here can be compared with Table 5 of Drasdo 2007,
% with the caveat that (e.g.) the superior meridian in the retina returned
% by this routine should be compared to the inferior meridian in the visual
% field presented in that table of Drasdo.



Rv0 = 0.011785;
Ro0 = 0.008333;

E20 = 20;
        if angle >= 0 && angle <= 90 
            E2v = interp1([0,90],[2.19,1.98],angle,interp);
         
        elseif angle > 90 && angle <= 180 
            E2v = interp1([90,180],[1.98,2.26],angle,interp);
          
        elseif angle > 180 && angle <= 270
            E2v = interp1([180,270],[2.26,1.5],angle,interp);
           
        elseif angle > 270 && angle <=360
            E2v = interp1([270,360],[1.5,2.19],angle,interp);
          
        end
Rve = Rv0.*(1+supportPosDeg./E2v);
Roe =Ro0.*(1+supportPosDeg./E20);

k=1+(1.004 - 0.007209.*r + 0.001694.*supportPosDeg.^2 - 0.00003765.*supportPosDeg.^3).^-2;

RFdensitySqDeg = ((1.12+0.0273.*r).*k)./(1.155.*(Rve.^2-Roe.^2));


RFdensitySqDeg=fliplr(RFdensitySqDeg);

    
end