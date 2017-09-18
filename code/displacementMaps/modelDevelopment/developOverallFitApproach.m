%% developOverallFitApproach
%
% This routine models retinal ganglion cell displacement.
% Our strategy is to begin with empirical measurements of cone and retinal
% ganglion cell densities obtained from each of the four cardinal meridians
% of the human retina. The data we use are from two papers published by
% Curcio and colleagues in 1990.
%
% We then engage in a modeling exercise to find a low-dimensional
% parameterization of the transformation of cone density into midget
% ganglion cell receptive field density (mRF) and of retinal ganglion cell
% density into midget retinal ganglion cell density (mRGC). In each case,
% the models do not make use of explicit information regarding the retinal
% position of the measurement to be transformed.
%
% We are then in a position to model mRGC and mRF density as a function of
% cone and RGC density, subject to a small number of parameters. Each of
% these functions in turn can be expressed as a cumulative mRGC and mRF
% count across the radial eccentricity of the retina.
%
% As described by Drasdo (2007), the cumulative counts of mRF and mRGCs
% will become equal at the eccentricity at which the RGCs are no longer
% displaced. Further, the mRF and mRGC cumulative counts should be
% equivalent beyond this point.
%
% We engage in a non-linear model fit of the parameters that transform
% cone density --> mRF density and RGC density --> mRGC density, subject
% to the constraint that the cumulative mRF density must be greater than
% the cumulative mRGC density within the displacement zone, and we minimize
% an error function defined as the difference in mRGC and mRF cumulative
% counts beyond the displacement zone. We set the displacement zone target
% to be 17 degrees eccentricity for all meridians, except for the nasal
% meridian within which the displacement point of 11 degrees is enforced by
% the presence of the optic nerve head.
%
% We find that a single set of parameters that governs the RGC --> mRGC
% transform is sufficient to model the data from all four meridians.
% Further, we find sets of parameters that vary only slightly between the
% meridians in the transform of cone density --> mRF density.
%



%% Housekeeping
close all
clear all
clc


%% Set up some variables

% Each of the meridians is defined by a polar angle value.
meridianNames = {'Nasal   ' 'Superior' 'Temporal' 'Inferior' };
meridianAngles = [0, 90, 180, 270 45];
meridianColors = {'g','b','r','k'};

% This is point in degrees at which displacement should become zero for
% each meridian
targetDisplacementPointDeg = [10 17 17 17];

% The calculations are performed across a regular sampling of eccentricity
% define a sample resolution. We note that the sample resolution must be
% sufficient fine so that the cumulative is an accurate estimate of the
% integral. Further, we find that our results depend in unpredictable ways
% on the particular maxModeledEccentricity selected. This latter value must
% be sufficiently outside the displacement zone so that there is a portion
% of the cumulative to match between the mRF and mRGC functions, but not so
% large as to venture into the periphery where our transform models are
% less accurate
sampleResolutionDegrees = 0.01;
maxModeledEccentricity = 30;
regularSupportPosDeg = 0:sampleResolutionDegrees:maxModeledEccentricity;


%% Derive parameters for the transformation of RGC density to mRGC density
[ rgcInitialTransformParams, figHandles(1) ] = developMidgetRGCFractionModel();


%% Derive parameters for the transformation of cone density to mRF density
[ rfInitialTransformParams, figHandles(2) ] = developMidgetRFFractionModel();

%% prepare for figures
figNames={'mRGCmodel','mRFmodel','displacement','mRGCwatsonCompare','mRFwatsonCompare'};
figOutpath='~/Desktop/displacementModelFigs';
if ~exist(figOutpath)
    mkdir(figOutpath)
end
figHandles(3) = figure;
figHandles(4) = figure;
figHandles(5) = figure;
for ff=1:length(figHandles)
    set(0, 'CurrentFigure', figHandles(ff))
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 1]);
end

%% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% mRF function
    % Build a function that returns mRF density over regular support.
    % We build the function using cone density, and subject to two fit
    % params
    
    % Obtain a spline fit to the empirical cone density data of Curcio 1990
    [coneDensityFit] = getSplineFitToConeDensity(meridianAngles(mm));
    
    % Create an anonymous function that returns mRF density as a function
    % cone density, with the transform defined by the first two fitParams
    mRFDensityOverRegularSupport = ...
        @(fitParams) transformConeToMidgetRFDensity(coneDensityFit(regularSupportPosDeg), ...
        'logitFitParams',fitParams(1:2))';
    % Define anonymous function for the cumulative sum of mRF density
    mRF_cumulative = @(fitParams) calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport(fitParams));
    
    
    %% mRGC function
    % Build a function that returns mRGC density over regular support.
    % We build the function using RGC density, and subject to the last
    % three fit params
    
    % Obtain a spline fit to the empirical RGC density data of Curcio 1990
    RGCDensityFit = getSplineFitToRGCDensity(meridianAngles(mm));
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
    set(0, 'CurrentFigure', figHandles(3))
    subplot(2,length(meridianAngles),mm);
    plot(regularSupportPosDeg(1:length(rgcDisplacementDeg)),rgcDisplacementDeg,'-r')
    ylim([-.5 3.0]);
    xlabel('eccentricity [deg]');
    ylabel('RGC displacement [deg]');
    title(meridianNames{mm});
    pbaspect([1 1 1]);
    
    % Plot the cumulative functions
    subplot(2,length(meridianAngles),mm+length(meridianAngles));
    plot(regularSupportPosDeg,mRGC_cumulative(fitParams(mm,:)),'-k')
    xlabel('eccentricity [deg]');
    ylabel('cells per sector');
    pbaspect([1 1 1]);
    hold on
    plot(regularSupportPosDeg,mRF_cumulative(fitParams(mm,:)),'-b')
    ylim([0 8e5]);
    if mm==1
        legend({'mRGC','mRF'},'Location','southeast');
    end
    hold off
    drawnow
    
    
    %% Plot the mRGC fraction
    set(0, 'CurrentFigure', figHandles(4))
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, RGCNativeSupportPosDeg ] = getCurcioRGCDensityByEccen( meridianAngles(mm) );
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
    [ ~, midgetFraction_ours ] = transformRGCToMidgetRGCDensity( RGCNativeSupportPosDeg', RGCDensitySqDeg', 'recipFitParams', fitParams(mm,3:5) );
    plot(RGCNativeSupportPosDeg,midgetFraction_ours,'-','Color',meridianColors{mm});
    hold on
    ylim([0 1]);
    xlim([0 40]);
    xlabel('eccentricity deg');
    ylabel('midget fraction');
    title('Our midget fraction');
    pbaspect([2 1 1]);
    drawnow


    %% Plots related to midget RF density
    set(0, 'CurrentFigure', figHandles(5))
    
    % load the empirical cone density measured by Curcio
    [coneDensitySqDeg, coneNativeSupportPosDeg] = getCurcioConeDensityByEccen(meridianAngles(mm));
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg));
    coneNativeSupportPosDeg = coneNativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    
    % calculate the mRF density using Watson equation 8 at the sites of
    % empirical cone measurement
    [ mRFDensitySqDeg_watson ] = calcWatsonMidgetRFDensityByEccen(coneNativeSupportPosDeg, meridianAngles(mm));
    
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
        'logitFitParams',fitParams(mm,1:2))';
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
    
    
end % loop over meridians

% Clean up the figures and save
set(0, 'CurrentFigure', figHandles(4))
subplot(1,2,2);
legend(meridianNames,'Location','southwest');

set(0, 'CurrentFigure', figHandles(5))
subplot(2,2,1);
legend(meridianNames,'Location','southwest');
subplot(2,2,2);
legend(meridianNames,'Location','southwest');

for ff = 1:length(figHandles)
    outFileName = fullfile(figOutpath,[figNames{ff} '.pdf']);
    saveas(figHandles(ff), outFileName, 'pdf');
end

% Dump the fitParams to the screen
outline='meridian\t\tslope\tinflect\tmRGC_a\tmRGC_b\tmRGC_c\n';
fprintf(outline);
for mm=1:length(meridianAngles)
    outline=[meridianNames{mm} '\t'];
    for pp=1:size(fitParams,2)
        outline=[outline '\t' num2str(fitParams(mm,pp))];
    end
    outline=[outline '\n'];
    fprintf(outline);
end

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

