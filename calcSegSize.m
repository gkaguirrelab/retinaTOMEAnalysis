function areaPerSeg= calcSegSize(radDeg,smpPerDeg,sectorAngle)
% 
% binWidthEcc =1/smpPerDeg;
% r_inner = 0:binWidthEcc:radDeg-binWidthEcc;
% r_outer = binWidthEcc:binWidthEcc:radDeg;
% 
% areaPerSeg =  ((pi.*r_outer.^2)-(pi.*r_inner.^2))./(360/sectorAngle);
% areaPerSeg = [0 areaPerSeg]
a=sectorAngle;
K=radDeg*smpPerDeg;
i_deg=1/smpPerDeg:1/smpPerDeg:radDeg+(1/smpPerDeg); 

%convert to mm 
ecc_mm = 0.268.*i_deg + 0.0003427.*(i_deg).*2 - (8.3309e-6).*(i_deg).^3;
alpha= 0.0752+5.846e-5*ecc_mm-1.064e-5*ecc_mm.^2+4.116e-8*ecc_mm.^3;

areaPerSeg = ((a*pi/360).*((2*ecc_mm.*(ecc_mm/2).^2)/K^2)).*(1./alpha);

end

