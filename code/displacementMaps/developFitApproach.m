
close all
clear all
clc


% define a sample resolution
sampleResolutionDegrees = 0.01;

% Each of the meridians is defined by a polar angle value.
meridianNames = {'Nasal' 'Superior' 'Temporal' 'Inferior'};
meridianAngles = [0, 90, 180, 270];

% These are the empirical displacement values found for each meridian. We
% want our function to produce these
targetDisplacementPointDeg = [12 17 17 17];

rgcCumProportionReferenceEccen = 20;

%% Step 0
% Examine the relationship between cone density and midget
% receptive field density. To do so, we load the Curcio cone density
% measurements for each meridian. At each of the sampled eccentricity
% values, we obtain the midgetRF density as provided by Eq 8 of Watson
% 2014.
%
% The goal is to define a function with few parameters that can relate cone
% density at any retinal location to midget receptive field density.
% We must find a transformation of the data such that the same fit function
% will well model the data for every meridian.
% A transformation that has this property is:
%      log10(log10(coneDensity) / mRFDensity)
% when plotted against log10(coneDensity), the relationship is well modeled
% by a 3 parameter reciprocal function. The temporal meridian has notably
% different parameter fits.

% Exclude from fitting those points from the far periphery
fitEccenLimit = 40; % degrees

% Set up some constants to deal with the log transform peculiarities
logMultiplierFix = 1e5;
logZeroFix = 1e-3;

% Set up a figure for this step
figHandle(1)=figure;

% Loop over the meridians
for mm = 1:length(meridianAngles)
    % load the empirical cone density measured by Curcio
    [supportPosDeg,coneDensitySqDeg] = curcioConeDensitySqDeg(meridianAngles(mm));
    % calculate the mRF density at these eccentricity locations using
    % Watson equation 8.
    [ midgetRFDensity_degSq ] = getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, meridianAngles(mm));
    
    % Set the foveal point to a non-zero value to permit the log transform
    supportPosDeg(1)=logZeroFix;
    
    % Remove nans and points beyond the modeled eccentricity bound
    isvalididx=find(~isnan(midgetRFDensity_degSq).*~isnan(coneDensitySqDeg) .* (supportPosDeg < fitEccenLimit) );
    supportPosDeg = supportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    midgetRFDensity_degSq = midgetRFDensity_degSq(isvalididx);
    
    % Define the transformed ratio function. Note that we multiply the
    % transformation by 1e5 to place the log transform in the non-zero
    % range
    ratioFunc = log10((log10(coneDensitySqDeg) ./ midgetRFDensity_degSq)*logMultiplierFix);
    x = log10(coneDensitySqDeg)';
    y = ratioFunc';
    
    % Set up a weight vector that forces the function to perfectly fit the
    % foveal value. This is needed so that the function returns the same
    % value at the fovea regardless of the meridian
    weights = ones(1,length(coneDensitySqDeg));
    weights(1) = 10000;
    
    % Define the reciprocal function and fit
    recipFunc = fittype('(1./(a+(b.*x)))+c','independent','x','dependent','y');
    recipFit = fit(x,y,recipFunc,'StartPoint',[-.1,.1,-2],'Weights',weights,'Lower',[-.2,0,-4],'Upper',[0,.2,-1] );
    
    % Save the fit params
    transformParams(mm).a=recipFit.a;
    transformParams(mm).b=recipFit.b;
    transformParams(mm).c=recipFit.c;
    
    % Plot the results
    subplot(2,length(meridianAngles),mm);
    plot(recipFit,x,y);
    legend off
    xlabel('log10 cone density per sq deg');
    ylabel('cone:mRF');
    title([meridianNames{mm} ' meridian']);
    
    subplot(2,length(meridianAngles),mm+length(meridianAngles));
    loglog(supportPosDeg,midgetRFDensity_degSq,'.k');
    hold on
    mRFDensityFit = transformConeToMidgetRFDensity(coneDensitySqDeg,[recipFit.a recipFit.b recipFit.c]);
    loglog(supportPosDeg,mRFDensityFit,'-b');
    xlabel('log10 eccentricity deg');
    ylabel('log 10 mRF density');
    
end
hold off

% prepare figures
figHandle(2) = figure;
figHandle(3) = figure;

% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% Step 1
    % Build a function that will return mRF denisty as a function of
    % eccentricity, taking as input the empirical cone density values and the
    % parameters that define the transform of cone density --> mRF density.
    
    % load the empirical cone density measured by Curcio
    [supportPosDeg,coneDensitySqDeg] = curcioConeDensitySqDeg(meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg)  );
    supportPosDeg = supportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    % Obtain a spline fit function for cone density
    coneDensityFit = fit(supportPosDeg',coneDensitySqDeg','smoothingspline','SmoothingParam', 1);
    % Switch to equal eccentricity sampling between 0 and 30 degrees
    regularSupportPosDeg = 0:sampleResolutionDegrees:30;
    % Create a display supports that handles the zero
    dispRegularSupportPosDeg=regularSupportPosDeg;
    dispRegularSupportPosDeg(1)=logZeroFix;
    dispSupportPosDeg=supportPosDeg;
    dispSupportPosDeg(1)=logZeroFix;
    % Create an anonymous function that returns mRF density as a function of
    % the parameters of the reciprocal transform
    mRFDensityOverRegularSupport = @(fitParams) transformConeToMidgetRFDensity(coneDensityFit(regularSupportPosDeg),fitParams(1:3))';
        
    
    %% Step 2
    % Build a function with the RGC density over the regular support
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, supportPosDeg ] = getCurcioRGCDensityByEccen( meridianAngles(mm) );
        % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    supportPosDeg = supportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    % Fit a spline to the RGC density data
    RGCDensityFit = fit(supportPosDeg,RGCDensitySqDeg,'smoothingspline', 'Exclude',find(isnan(RGCDensitySqDeg)),'SmoothingParam', 1);
    % Create an anonymous function that returns mRGC density as a function of
    % RGC density and the transformParams
    mRGCDensityOverRegularSupport = @(fitParams) transformRGCToMidgetRGCDensity(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg)','recipFitParams',fitParams(4:6));

    
    %% Step 3
    % Our goal is to have the cumulative mRF and mRGC functions have minimally
    % different values past the displacement point, and for the cumulative mRF
    % to have a greater value than the mRGC function prior to the displacement point
    
    % Set up an anonymous function and a variable that will hold the mRF and
    % mRGC cumulative sums, respectively
    mRF_ringcount = @(fitParams) calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport(fitParams));
    mRGC_ringcount = @(fitParams)calcCumulative(regularSupportPosDeg, mRGCDensityOverRegularSupport(fitParams));
    % Create a non-linear constraint that tests if the RF values are greater
    % than the RGC values at eccentricities less than the displacement point
    nonlinconst = @(fitParams) testRFGreaterThanRGC(regularSupportPosDeg, mRF_ringcount(fitParams), mRGC_ringcount(fitParams), targetDisplacementPointDeg(mm));
    % Define an error function
    errorFunc = @(fitParams) errorMatchingRFandRGC(regularSupportPosDeg, mRF_ringcount(fitParams), mRGC_ringcount(fitParams), targetDisplacementPointDeg(mm));
    % set some bounds
    lb = [-0.2 0 -4 2.4026 -8.0877 -0.0139];
    ub = [0.0 .2 -1 2.4026 -8.0877 -0.0139];
    x0 = [-0.1421    0.1370   -1.1503 2.4026 -8.0877 -0.0139];
    % Run the fitter
    fitParams(mm,:) = fmincon(errorFunc,x0,[],[],[],[],lb,ub,nonlinconst);
    % Plot the ringcounts
    set(0, 'CurrentFigure', figHandle(2))
    subplot(2,length(meridianAngles),mm+length(meridianAngles));
    plot(regularSupportPosDeg,mRGC_ringcount(fitParams(mm,:)),'-k')
    xlabel('eccentricity [deg]');
    ylabel('cells per sector');
    hold on
    plot(regularSupportPosDeg,mRF_ringcount(fitParams(mm,:)),'-b')
    ylim([0 8e5]);
    title('mRGC black, mRF blue');
    hold off
    
    
    %% Step 4
    % Calculate displacement and plot
    rgcDisplacementDeg=calcDisplacement(regularSupportPosDeg, mRGC_ringcount(fitParams(mm,:)), mRF_ringcount(fitParams(mm,:)));
    set(0, 'CurrentFigure', figHandle(2))
    subplot(2,length(meridianAngles),mm);
    plot(regularSupportPosDeg(1:length(rgcDisplacementDeg)),rgcDisplacementDeg,'-r')
    ylim([-.5 3.0]);
    xlabel('eccentricity [deg]');
    ylabel('radial displacement of RGCs [deg]');
    title([meridianNames{mm} ' meridian']);
    
    
    %% Step 5
    % Compare the midgetRF : cone ratio for each meridian between the
    % Watson model and what we find is needed for good dispacement
    
    % load the empirical cone density measured by Curcio
    [supportPosDeg,coneDensitySqDeg] = curcioConeDensitySqDeg(meridianAngles(mm));
    % calculate the mRF density at these eccentricity locations using
    % Watson equation 7.
    [ midgetRFDensity_degSq ] = getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, meridianAngles(mm));

    % Get the log-friendly x axis values ready
    dispSupportPosDeg=supportPosDeg;
    dispSupportPosDeg(1)=logZeroFix;

    % Plot this ratio
    set(0, 'CurrentFigure', figHandle(3))
    subplot(1,length(meridianAngles),mm);
        
    loglog(dispSupportPosDeg,midgetRFDensity_degSq./coneDensitySqDeg,'.k');
    ylim([1e-4 1e1]);
    hold on
    xlabel('log10 eccentricity');
    ylabel('log10 cone:mRF');
    title([meridianNames{mm} ' - Watson (.), Adjusted (-)']);

    % Get the midget receptive field density implied by the found transform
    % parameters
    mRFDensityDegSq_adjusted=transformConeToMidgetRFDensity(coneDensityFit(supportPosDeg)',fitParams(mm,:));

    % Plot this ratio
    loglog(dispSupportPosDeg,mRFDensityDegSq_adjusted./coneDensitySqDeg,'-r');
    hold off
    
    

    
end % loop over meridians




% function mRGCDensity = transformRGCToMidgetRGCDensity(regularSupportPosDeg,rgcDensitySqDeg,transformParams)
% 
% % Distribute the transformParams
% rm = transformParams(4);
% 
% % Calculate the midget fraction.
% f0 = 0.8928; % One of Watson's constants
% 
% midgetFraction = calcWatsonMidgetFractionByEccen(regularSupportPosDeg,f0,rm);
% 
% mRGCDensity = rgcDensitySqDeg .* midgetFraction;
% 
% end


function mRFDensity = transformConeToMidgetRFDensity(coneDensitySqDeg,transformParams)

% Distribute the transformParams
a = transformParams(1);
b = transformParams(2);
c = transformParams(3);

mRFDensity = log10(coneDensitySqDeg) ./ (10.^( ( 1./(a+(b.*log10(coneDensitySqDeg))) )+c-5));

end




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

