
% clear out the workspace
close all
clear all
clc

% make the displacement map using the default params
fprintf('*** makeDisplacementMap\n');

[ displacementMapDeg, fitParams, meridianAngles, rgcDisplacementEachMeridian, mRGC_cumulativeEachMeridian, mRF_cumulativeEachMeridian ] = ...
    makeDisplacementMap('verbose',true);

fprintf('*** done\n');

% Plot the displacement map
figure
maxDisplacementDeg = max(displacementMapDeg(:));
climVals = [0,ceil(maxDisplacementDeg)];
imagesc(displacementMapDeg, climVals);
c = colorbar;
c.Label.String='Radial RGC displacement [deg]';
xlabel('Position [deg] temporal --> nasal');
ylabel('Position [deg] inferior --> superior');
numTicks=length(xticks);
k=linspace(-1*30,30,numTicks+1);
xticklabels(string(k(2:end)))
yticklabels(string(k(2:end)))

% Plot the cumulatives and displacements across meridians
figure
regularSupportPosDeg = 0:0.01:30;
for mm = 1:length(meridianAngles)
    
    % plot the displacement
    subplot(length(meridianAngles),2,mm*2);
    plot(regularSupportPosDeg,rgcDisplacementEachMeridian(mm,:),'-r')
    axis off;
    ylim([-.5 3.0]);
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('RGC displacement [deg]');
    end
    
    % Plot the cumulative functions
    subplot(length(meridianAngles),2,mm*2-1);
    plot(regularSupportPosDeg,mRGC_cumulativeEachMeridian(mm,:),'-k')
    axis off;
    if mm == length(meridianAngles)
        axis on;
        xlabel('eccentricity [deg]');
        ylabel('cells per sector');
    end
    hold on
    plot(regularSupportPosDeg,mRF_cumulativeEachMeridian(mm,:),'-b')
    ylim([0 8e5]);
    hold off
    drawnow
end

% Show the cone --> mRF model
developMidgetRFFractionModel( 'makePlots', true );

% Show the RGC --> mRGC model
developMidgetRGCFractionModel( 'makePlots', true );

% Plot the spline cone density fits
figure
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};
subplot(1,2,1)
for mm=1:4
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(cardinalMeridianAngles(mm));
    dispConeNativeSupportPosDeg=coneNativeSupportPosDeg;
    dispConeNativeSupportPosDeg(1)=1e-2;
    [coneDensityFit] = getSplineFitToConeDensity(cardinalMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
end
legend(num2str(cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [coneDensityFit] = getSplineFitToConeDensity(interpolarMeridianAngles(mm));
    loglog(dispConeNativeSupportPosDeg,coneDensityFit(coneNativeSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 Cone density [counts / deg2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')


% Plot the spline RGC density fits
figure
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};
subplot(1,2,1)
for mm=1:4
    [RGCDensitySqDeg, RGCNativeSupportPosDeg] = getCurcioRGCDensityByEccen(cardinalMeridianAngles(mm));
    dispRGCNativeSupportPosDeg=RGCNativeSupportPosDeg;
    dispRGCNativeSupportPosDeg(1)=1e-2;
    [RGCDensityFit] = getSplineFitToRGCDensity(cardinalMeridianAngles(mm));
    loglog(dispRGCNativeSupportPosDeg,RGCDensitySqDeg,'x','Color',meridianColors{mm});
    xlim([1e-2,70]);
    hold on
    regularSupportPosDeg=1e-2:0.01:70;
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
end
legend(num2str(cardinalMeridianAngles),'Location','southwest')

subplot(1,2,2)
interpolarMeridianAngles=[45 135 225 315];
for mm=1:4
    [RGCDensityFit] = getSplineFitToRGCDensity(interpolarMeridianAngles(mm));
    loglog(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg),'-','Color',meridianColors{mm});
    xlim([1e-2,70]);
    xlabel('log10 Eccentricity [deg]');
    ylabel('log10 RGC density [counts / deg2]');
    hold on
end
legend(num2str(interpolarMeridianAngles),'Location','southwest')


%% Plot the mRGC fraction for the cardinal meridians
figure
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};

for mm = 1:4
    
    meridianIdx = find(meridianAngles==cardinalMeridianAngles(mm),1);
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = getCurcioRGCDensityByEccen( cardinalMeridianAngles(mm) );
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    RGCNativeSupportPosDeg = RGCNativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    
    % Plot Watson's midget fraction
    subplot(1,2,1);
    f0 = 0.8928; rm = 41.03; % Watson's values
    midgetFraction_watson = calcWatsonMidgetFractionByEccen(RGCNativeSupportPosDeg,f0,rm);
    plot(RGCNativeSupportPosDeg,midgetFraction_watson,'-k');
    hold on
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    ylim([0 1]);
    xlim([0 40]);
    title('Watson''s midget fraction (from Drasdo)');
    pbaspect([2 1 1]);
    
    % Plot our midget fraction
    subplot(1,2,2);
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensity( RGCNativeSupportPosDeg', RGCDensitySqDeg', 'recipFitParams', fitParams(meridianIdx,3:5) );
    plot(RGCNativeSupportPosDeg,midgetFraction_ours,'-','Color',meridianColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    title('Our midget fraction');
    pbaspect([2 1 1]);
    drawnow
end


%% Plots related to midget RF density

figure
cardinalMeridianAngles=[0 90 180 270];
meridianColors={'k','r','g','b'};

for mm = 1:4
    
    meridianIdx = find(meridianAngles==cardinalMeridianAngles(mm),1);
    
    % load the empirical cone density measured by Curcio
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(cardinalMeridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    
    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDeg_watson ] = calcWatsonMidgetRFDensityByEccen(coneNativeSupportPosDeg, cardinalMeridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Watson''s mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDeg_ours = transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg), ...
        'logitFitParams',fitParams(meridianIdx,1:2))';
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Our mRF density by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e1]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    ylim([1e-4 1e2]);
    title('Watson''s mRF:cone ratio by eccentricity');
    pbaspect([2 1 1]);
    hold on
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    mRFDensitySqDeg_ours = ...
        transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg)','logitFitParams',fitParams(mm,1:2));
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e2]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    title('Our mRF:cone ratio by eccentricity');
    hold on
    pbaspect([2 1 1]);
    drawnow
    
end