radDeg = 20;
smpPerDeg=2;
sectorAngle = 6;

%% Generate the RGCf Density
retData = densityRGC(radDeg,smpPerDeg,'OFF');

%% Generate the RGCf Density
retData(:,:,4) = densityRf(radDeg,smpPerDeg,'OFF');



rotDeg = 0; %nasal
displace0 = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg);


rotDeg = 90; % inf
displace90 = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg);



rotDeg = 180; % temp 
displace180 = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg);

rotDeg = 270; % sup
displace270 = calcDisp(retData,radDeg,smpPerDeg,sectorAngle,rotDeg);



figure 
hold on 
plot(0:1/smpPerDeg:radDeg,displace0,'g'); %nasal
plot(0:1/smpPerDeg:radDeg,displace90,'k'); %inf
plot(0:1/smpPerDeg:radDeg,displace180,'r')%temp
plot(0:1/smpPerDeg:radDeg,displace270,'b'); %sup
legend('nasal','inferior','temporal','superior');
