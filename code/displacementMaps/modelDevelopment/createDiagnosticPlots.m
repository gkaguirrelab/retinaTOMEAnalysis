close all

% Plot the cone density fits
figure
meridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(meridianAngles(mm));
    dispConeNativeSupportPosDeg=coneNativeSupportPosDeg;
    dispConeNativeSupportPosDeg(1)=1e-2;
    [coneDensityFit] = getSplineFitToConeDensity(meridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
end
legend(num2str(meridianAngles),'Location','southwest')

subplot(1,2,2)
meridianAngles=[45 135 225 315];
for mm=1:4
    [coneDensityFit] = getSplineFitToConeDensity(meridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
    hold on
end
legend(num2str(meridianAngles),'Location','southwest')


% Plot the RGC density fits
figure
meridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};
subplot(1,2,1)
for mm=1:4
    [RGCDensitySqDeg, RGCNativeSupportPosDeg] = getCurcioRGCDensityByEccen(meridianAngles(mm));
    dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
    dispRGCNativeSupportPosDeg(1)=1e-2;
    [RGCDensityFit] = getSplineFitToRGCDensity(meridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    regularSupportPosDeg=1e-2:0.01:70;
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
end
legend(num2str(meridianAngles),'Location','southwest')

subplot(1,2,2)
meridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensity(meridianAngles(mm));
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
    hold on
end
legend(num2str(meridianAngles),'Location','southwest')






