radDeg=20;
smpPerDeg = 2;

midPt = round((2*(radDeg*smpPerDeg))/2)+1;
figure
hold on 
xRange = 0:1/smpPerDeg:radDeg;
xRangeBck = radDeg:-1/smpPerDeg:0;
plot(xRangeBck,meridian(midPt,1:midPt,4),'r')%temp

plot(xRangeBck,meridian(1:midPt,midPt,4),'b'); %sup

plot(xRange,meridian(midPt,midPt:end,4),'g'); %nasal

plot(xRange,meridian(midPt:end,midPt,4),'k'); %inf

set(gca,'xscale','log')
set(gca,'yscale','log')
