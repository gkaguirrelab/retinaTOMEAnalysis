
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
targetDisplacementPointDeg = [12 15 15 15];
% The x-axis value we will assign "zero" for log transforms
logZeroPlotFix = 1e-5;


%% Derive parameters for the transformation of RGC density to mRGC density
[ rgcInitialTransformParams, figureHandle(1) ] = developMidgetRGCFractionModel();


%% Derive parameters for the transformation of cone density to mRF density
[ rfInitialTransformParams, figureHandle(2) ] = developMidgetRFFractionModel();


%% prepare figure handles
figHandle(3) = figure;
figHandle(4) = figure;


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

    
    %% Plot the results
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
    
    % Compare the midgetRF : cone ratio for each meridian between the
    % Watson model and what we find is needed for good dispacement
    
    % calculate the mRF density at locations for which we have empirical
    % cone measurements using Watson equation 8.
    [ mRFDensitySqDeg_watson ] = getWatsonMidgetReceptiveFieldDensityByEccen(coneNativeSupportPosDeg, meridianAngles(mm));

    % Get the log-friendly x axis values ready
    dispSupportPosDeg=coneNativeSupportPosDeg;
    dispSupportPosDeg(1)=logZeroPlotFix;

    % Plot this ratio
    set(0, 'CurrentFigure', figHandle(4))
    subplot(1,2,1);
    loglog(dispSupportPosDeg,mRFDensitySqDeg_watson./coneDensitySqDeg,'-','Color',meridianColors{mm});
    ylim([1e-4 1e1]);
    xlabel('log10 eccentricity');
    ylabel('log10 mRF:cone');
    ylim([1e-4 1e2]);
    title('mRF:cone ratio by eccentricity from Watson');
    hold on
    
    % Get the midget receptive field density implied by the found transform
    % parameters
    mRFDensitySqDeg_ours=transformConeToMidgetRFDensity(coneDensityFit(coneNativeSupportPosDeg)','logitFitParams',fitParams(mm,1:2));

    % Plot this ratio
    subplot(1,2,2);
    loglog(dispSupportPosDeg,mRFDensitySqDeg_ours./coneDensitySqDeg,'-','Color',meridianColors{mm});
    ylim([1e-4 1e2]);
    hold on
    drawnow
    
end % loop over meridians

% Clean up some figure legends
set(0, 'CurrentFigure', figHandle(4))
subplot(1,2,1);
legend(meridianNames,'Location','southwest');
subplot(1,2,2);
legend(meridianNames,'Location','southwest');


function [c,ceq] = testRFGreaterThanRGC(regularSupportPosDeg, countsPerRingRF, countsPerRingRGC, displacementPointDeg)

% If there are any cumulative RGC values that are greater than the RF
% values at eccentricities less than the displacementPoint, then this
% violates the nonlinear fit constraint

withinRangeIdx = find(regularSupportPosDeg < displacementPointDeg);
c = sum(countsPerRingRGC(withinRangeIdx) > countsPerRingRF(withinRangeIdx));
ceq = []; % unused
end

function error = errorMatchingRFandRGC(regularSupportPosDeg, countsPerRingRF, countsPerRingRGC, displacementPointDeg)

% The error is calculated as the difference between the cumulative RF and
% RGC counts at retinal positions past the point where displacement should
% have ended

withinRangeIdx = find(regularSupportPosDeg > displacementPointDeg);
error = sqrt(sum((countsPerRingRGC(withinRangeIdx) - countsPerRingRF(withinRangeIdx)).^2));
end


function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

tmp = diff(regularSupportPosDeg);
sampleResolutionDegrees = tmp(1);

% Measure the displacement (in degrees)
displaceInSamples=arrayfun(@(x) find(countPerRingRGC>x,1), countPerRingRF(2:end),'UniformOutput',false);
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
displaceInDeg = (displaceInSamples - (1:1:length(displaceInSamples))) * sampleResolutionDegrees;
displaceInDeg(find(displaceInDeg <= 0,1):end)=0;

end

