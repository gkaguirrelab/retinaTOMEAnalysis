function validatePlots(RGCdensity,RFdensity,countRF,countRGC,radDeg,smpPerDeg)
mdPt = round(size(RGCdensity,1)/2);
figure; 

subplot(2,2,1)% plot the HM and VM cross section of the RGC denstiy image
hold on
xSmps = -radDeg:1/smpPerDeg:radDeg;
plot(xSmps,RGCdensity(mdPt,:),'r')
plot(xSmps,RGCdensity(:,mdPt),'b')
legend('Temporal-Nasal','Superior-Inferior')
xlabel('Eccentricity (deg)')
ylabel('RGC Denstiy (deg^{-2}')
title('RGC Denstiy by Eccentricity')
subplot(2,2,2)% plot the HM and VM cross section of the RF denstiy image
hold on
plot(xSmps,RFdensity(mdPt,:),'r')
plot(xSmps,RFdensity(:,mdPt),'b')
legend('Temporal-Nasal','Superior-Inferior')
xlabel('Eccentricity (deg)')
ylabel('RF Denstiy (deg^{-2}')
title('RF Denstiy by Eccentricity')
subplot(2,2,3) % plot Nasal cross section of the RGC & RF denstiy adjusted by area of subsector 
hold on
xSmps = 1/smpPerDeg:1/smpPerDeg:radDeg+1/smpPerDeg;
plot(xSmps,countRF,'r')
plot(xSmps,countRGC,'b')
legend('Receptive Feild','Ganglion Cell')
xlabel('Eccentricity (deg)')
ylabel('Cell Count in a Subsector')
title('RGC & RF Cell Count in a Subsector')
subplot(2,2,4) % plot Nasal cross section of the RGC & RF cumulative cell count from the subsectors 
hold on
xSmps = 1/smpPerDeg:1/smpPerDeg:radDeg+1/smpPerDeg;
plot(xSmps,cumsum(countRF),'r')
plot(xSmps,cumsum(countRGC),'b')
legend('Receptive Feild','Ganglion Cell')
xlabel('Eccentricity (deg)')
ylabel('Cumulative Cell Count')
title('Cumulative Cell Count Across Subsectors')
end