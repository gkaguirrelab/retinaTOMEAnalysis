radMM = 5;
smpPerMM=4;
sectorAngle = 6;

%% Generate the RGCf Density
retData = rgcDensityMM(radMM,smpPerMM);

%% Generate the RGCf Density
retData(:,:,4) = rfDensityMM(radMM,smpPerMM);



rotDeg = 0; %nasal
displace0 = calcDisp(retData,radMM,smpPerMM,sectorAngle,rotDeg);


rotDeg = 90; % inf
displace90 = calcDisp(retData,radMM,smpPerMM,sectorAngle,rotDeg);



rotDeg = 180; % temp 
displace180 = calcDisp(retData,radMM,smpPerMM,sectorAngle,rotDeg);

rotDeg = 270; % sup
displace270 = calcDisp(retData,radMM,smpPerMM,sectorAngle,rotDeg);



figure 
hold on 
plot(0:1/smpPerMM:radMM,displace0,'g'); %nasal
plot(0:1/smpPerMM:radMM,displace90,'k'); %inf
plot(0:1/smpPerMM:radMM,displace180,'r')%temp
plot(0:1/smpPerMM:radMM,displace270,'b'); %sup
legend('nasal','inferior','temporal','superior');
