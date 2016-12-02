function areaPerSeg= calcSegSize(radMM,smpPerMM,sectorAngle)
% % 
% binWidthMM =1/smpPerMM;
% r_inner = binWidthMM:binWidthMM:radMM+binWidthMM;
% r_outer = binWidthMM+binWidthMM:binWidthMM:radMM+2*binWidthMM;
% 
% areaPerSeg =  ((pi.*r_outer.^2)-(pi.*r_inner.^2))./(360/sectorAngle);

a=sectorAngle;
K=radMM*smpPerMM;
ecc_mm=1/smpPerMM:1/smpPerMM:radMM+(1/smpPerMM); 

areaPerSeg = ((a*pi/360).*((2*ecc_mm.*(ecc_mm/2).^2)/K^2));

end

