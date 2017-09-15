
%% Housekeeping
close all
clear all
clc


%% Setup some variables
% The calculations are performed across a regular sampling of eccentricity
% define a sample resolution
sampleResolutionDegrees = 0.01;
maxModeledEccentricity = 30;
regularSupportPosDeg = 0:sampleResolutionDegrees:maxModeledEccentricity;
% Each of the meridians is defined by a polar angle value.
meridianNames = {'Nasal' 'Superior' 'Temporal' 'Inferior'};
meridianAngles = [0, 90, 180, 270];
meridianColors = {'g','b','r','k'};
% This is point in degrees at which displacement should become zero for
% each meridian
targetDisplacementPointDeg = [11 17 17 17];


%% Derive parameters for the transformation of RGC density to mRGC density
[ rgcInitialTransformParams, figureHandle(1) ] = developMidgetRGCFractionModel();


%% Derive parameters for the transformation of cone density to mRF density
[ rfInitialTransformParams, figureHandle(2) ] = developMidgetRFFractionModel();


%% prepare figure handles
figHandle(3) = figure;
figHandle(4) = figure;
figHandle(5) = figure;

%% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% mRF function
    % Build a function that returns mRF density over regular support.
    % We build the function using cone density.
    
    % load the empirical cone density measured by Curcio
    [coneNativeSupportPosDeg,coneDensitySqDeg] = getCurcioConeDensitySqDeg(meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    % Obtain a spline fit function for cone density
    coneDensityFit = fit(coneNativeSupportPosDeg',coneDensitySqDeg','smoothingspline','SmoothingParam', 1);
    % Create an anonymous function that returns mRF density as a function
    % cone density, with the transform defined by the first two fitParams
    mRFDensityOverRegularSupport = ...
        @(fitParams) transformConeToMidgetRFDensity(coneDensityFit(regularSupportPosDeg), ...
        'logitFitParams',fitParams(1:2))';
    % Define anonymous function for the cumulative sum of mRF density
    mRF_cumulative = @(fitParams) calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport(fitParams));
    
    
    %% mRGC function
    % Build a function that returns mRGC density over regular support.
    % We build the function using RGC density.
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = getCurcioRGCDensityByEccen( meridianAngles(mm) );
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    RGCNativeSupportPosDeg = RGCNativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    % Fit a spline to the RGC density data
    RGCDensityFit = fit(RGCNativeSupportPosDeg,RGCDensitySqDeg,'smoothingspline', 'Exclude',find(isnan(RGCDensitySqDeg)),'SmoothingParam', 1);
    % Create an anonymous function that returns mRGC density as a function of
    % RGC density, with the transform defined by the last three fitParams
    mRGCDensityOverRegularSupport = ...
        @(fitParams) transformRGCToMidgetRGCDensity(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg)',...
        'recipFitParams',fitParams(3:5));
    % Define anonymous function for the cumulative sum of mRF density
    mRGC_cumulative = @(fitParams) calcCumulative(regularSupportPosDeg, mRGCDensityOverRegularSupport(fitParams));
    
    
    %% Non-linear constraint and error functions
    % Our goal is to have the cumulative mRF and mRGC functions have minimally
    % different values past the displacement point, and for the cumulative mRF
    % to have a greater value than the mRGC function prior to the displacement point
    
    % Create a non-linear constraint that tests if the RF cumulative values
    % are greater than the RGC cumulative values at eccentricities less
    % than the displacement point
    nonlinconst = @(fitParams) testRFGreaterThanRGC(regularSupportPosDeg, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), targetDisplacementPointDeg(mm));
    
    % Define an error function
    errorFunc = @(fitParams) errorMatchingRFandRGC(regularSupportPosDeg, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), targetDisplacementPointDeg(mm));
    
    
    %% Perform the fit
    % We will search over the mRF transform parameters and lock the mRGC
    % transform parameters. Set upper and lower bounds on the mRF params
    % to be 1.25x the median values found across meridians (with a bit of
    % sign exponent trickery to handle the direction of negative params)
    lb = [rfInitialTransformParams./(1.25.^sign(rfInitialTransformParams)) rgcInitialTransformParams];
    ub = [rfInitialTransformParams.*(1.25.^sign(rfInitialTransformParams)) rgcInitialTransformParams];
    x0 = [rfInitialTransformParams rgcInitialTransformParams];
    
    % Set up the options
    options = optimoptions('fmincon','Display','none');
    
    % Fit that sucker
    fitParams(mm,:) = fmincon(errorFunc,x0,[],[],[],[],lb,ub,nonlinconst,options);
    
    
    %% Calculate displacement
    rgcDisplacementDeg=calcDisplacement(regularSupportPosDeg, mRGC_cumulative(fitParams(mm,:)), mRF_cumulative(fitParams(mm,:)));
    
    
    %% Plot the displacement and cumulative functions
    % plot the displacement
    set(0, 'CurrentFigure', figHandle(3))
    subplot(2,length(meridianAngles),mm);
    plot(regularSupportPosDeg(1:length(rgcDisplacementDeg)),rgcDisplacementDeg,'-r')
    ylim([-.5 3.0]);
    xlabel('eccentricity [deg]');
    ylabel('radial displacement of RGCs [deg]');
    title([meridianNames{mm} ' meridian']);
    
    % Plot the cumulative functions
    subplot(2,length(meridianAngles),mm+length(meridianAngles));
    plot(regularSupportPosDeg,mRGC_cumulative(fitParams(mm,:)),'-k')
    xlabel('eccentricity [deg]');
    ylabel('cells per sector');
    hold on
    plot(regularSupportPosDeg,mRF_cumulative(fitParams(mm,:)),'-b')
    ylim([0 8e5]);
    title('mRGC black, mRF blue');
    hold off
    drawnow
    
    
    %% Plots related to midget RF density
    set(0, 'CurrentFigure', figHandle(4))
    
    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDeg_watson ] = getWatsonMidgetReceptiveFieldDensityByEccen(coneNativeSupportPosDeg, meridianAngles(mm));
    
    % Plot the mRF density by eccentricity for Watson
    subplot(2,2,1);
    
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Watson''s mRF density by eccentricity');
    hold on
    
    % Plot the mRF density by eccentricity from our functions
    subplot(2,2,2);
    mRFDensitySqDeg_ours = transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg), ...
        'logitFitParams',fitParams(mm,1:2))';
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end),'-','Color',meridianColors{mm});
    ylim([1e0 1e5]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF density / deg2');
    title('Our mRF density by eccentricity');
    hold on
    
    % Plot the mRF : cone ratio for Watson
    subplot(2,2,3);
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_watson(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e1]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    ylim([1e-4 1e2]);
    title('Watom''s mRF:cone ratio by eccentricity');
    hold on
    
    % Plot the mRF : cone ratio for us
    subplot(2,2,4);
    mRFDensitySqDeg_ours = ...
        transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg)','logitFitParams',fitParams(mm,1:2));
    loglog(coneNativeSupportPosDeg(2:end),mRFDensitySqDeg_ours(2:end)./coneDensitySqDeg(2:end),'-','Color',meridianColors{mm});
    ylim([1e-4 1e2]);
    title('Our mRF:cone ratio by eccentricity');
    hold on
    drawnow
    
    
    %% Plot the mRGC fraction
    set(0, 'CurrentFigure', figHandle(5))
    
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
    
    % Plot our midget fraction
    subplot(1,2,2);
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensity( RGCNativeSupportPosDeg', RGCDensitySqDeg', 'recipFitParams', fitParams(mm,3:5) );
    plot(RGCNativeSupportPosDeg,midgetFraction_ours,'-','Color',meridianColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    title('Our midget fraction');
    drawnow
    
end % loop over meridians

% Clean up some figure legends
set(0, 'CurrentFigure', figHandle(4))
subplot(2,2,1);
legend(meridianNames,'Location','southwest');
subplot(2,2,2);
legend(meridianNames,'Location','southwest');

set(0, 'CurrentFigure', figHandle(5))
subplot(1,2,2);
legend(meridianNames,'Location','southwest');




%% LOCAL FUNCTIONS

function [c,ceq] = testRFGreaterThanRGC(regularSupportPosDeg, countsPerRingRF, countsPerRingRGC, displacementPointDeg)

% If there are any cumulative RGC values that are greater than the RF
% values at eccentricities less than the displacementPoint, then this
% violates the nonlinear fit constraint

withinRangeIdx = find(regularSupportPosDeg < displacementPointDeg);
c = sum(countsPerRingRGC(withinRangeIdx) > countsPerRingRF(withinRangeIdx));
ceq = []; % unused
end

function error = errorMatchingRFandRGC(regularSupportPosDeg, countsPerRingRF, countsPerRingRGC, displacementPointDeg)

% The error is calculated as the SSQ of the absolute difference between the
% cumulative RF and RGC counts at retinal positions past the point where
% displacement should have ended

withinRangeIdx = find(regularSupportPosDeg > displacementPointDeg);
error = sqrt(sum((countsPerRingRGC(withinRangeIdx) - countsPerRingRF(withinRangeIdx)).^2));
end

function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDeg);
sampleResolutionDegrees = tmp(1);

% Measure the displacement (in degrees)
displaceInSamples=arrayfun(@(x) find(countPerRingRGC>x,1), countPerRingRF(2:end),'UniformOutput',false);
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
displaceInDeg = (displaceInSamples - (1:1:length(displaceInSamples))) * sampleResolutionDegrees;
displaceInDeg(find(displaceInDeg <= 0,1):end)=0;

end

