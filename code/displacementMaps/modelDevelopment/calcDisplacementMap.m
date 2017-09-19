function [ displacementMapDeg, fitParams, meridianAngles ] = calcDisplacementMap( varargin )
% calcDisplacementMap( varargin )
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

% OUTPUT
%   displacementMapDeg -
%   fitParams - The values of the five parameters that adjust the transform
%      of cone --> mRF and RGC --> mRGC, provided for each meridian
%   meridianAngles - a vector of polar angle values (in degrees) for which
%       the displacement values were calculated
%
% OPTIONS
%   sampleResolutionDegrees, maxModeledEccentricity - The calculations are
%       performed across a regular sampling of eccentricity. These params
%       deine sample resolution and the max modeled eccentricity. We
%       note that the sample resolution must be sufficient fine so that the
%       cumulative is an accurate estimate of the integral. Further, we
%       find that our results depend in unpredictable ways on the
%       particular maxModeledEccentricity selected. This latter value must
%       be sufficiently outside the displacement zone so that there is a
%       portion of the cumulative to match between the mRF and mRGC
%       functions, but not so large as to venture into the periphery where
%       our transform models areless accurate
%   targetDisplacementPointDeg - This is point in degrees at which
%       displacement should become zero for each cadinal meridian
%   meridianAngleResolutionDeg - The resolution across polar angle for
%       which displacements are calculated.
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% Optional anaysis params
p.addParameter('sampleResolutionDegrees',0.01,@isnumeric);
p.addParameter('maxModeledEccentricity',30,@isnumeric);
p.addParameter('targetDisplacementPointDeg',[11 17 17 17 17 17 17 17],@isnumeric);
p.addParameter('meridianAngleResolutionDeg',45,@isnumeric);


% Optional display params
p.addParameter('verbose',true,@islogical);
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(varargin{:})


%% Setup
% Prepare the regular eccentricity support base
regularSupportPosDeg = ...
    0:p.Results.sampleResolutionDegrees:p.Results.maxModeledEccentricity;

% Prepare the set of meridian angles for which we will calculate
% displacement
meridianAngles = 0:p.Results.meridianAngleResolutionDeg:(360-p.Results.meridianAngleResolutionDeg);

% Create a figure to hold the displacement profiles
if p.Results.makePlots
    figure
end


%% Derive parameters for the transformation of RGC density to mRGC density
[ rgcInitialTransformParams ] = developMidgetRGCFractionModel('makePlots',false);


%% Derive parameters for the transformation of cone density to mRF density
[ rfInitialTransformParams ] = developMidgetRFFractionModel('makePlots',false);


%% Loop over the meridians
for mm = 1:length(meridianAngles)
    
    %% mRF function
    % We build a function that returns the cumulative mRF density, subject
    % to two fit parameters. This function is based upon a model of cone
    % density.
    
    % Obtain a spline fit to the empirical cone density data of Curcio 1990
    [coneDensityFit] = getSplineFitToConeDensity(meridianAngles(mm));
    % Create an anonymous function that returns mRF density as a function
    % of cone density, with the transform defined by the first two fitParams
    mRFDensityOverRegularSupport = ...
        @(fitParams) transformConeToMidgetRFDensity(coneDensityFit(regularSupportPosDeg), ...
        'logitFitParams',fitParams(1:2))';
    % Define anonymous function for the cumulative sum of mRF density
    mRF_cumulative = @(fitParams) calcCumulative(regularSupportPosDeg, mRFDensityOverRegularSupport(fitParams));
    
    
    %% mRGC function
    % We build a function that returns the cumulative mRGC density, subject
    % to three fit parameters. This function is based upon a model of RGC
    % density.
    
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
    % Create a non-linear constraint that tests if the RF cumulative values
    % are greater than the RGC cumulative values at eccentricities less
    % than the displacement point
    nonlinconst = @(fitParams) testRFGreaterThanRGC(regularSupportPosDeg, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), p.Results.targetDisplacementPointDeg(mm));
    
    % The error function acts to minimize the diffrence between the
    % mRF and mRGC cumulative functions past the displacement point
    errorFunc = @(fitParams) errorMatchingRFandRGC(regularSupportPosDeg, mRF_cumulative(fitParams), mRGC_cumulative(fitParams), p.Results.targetDisplacementPointDeg(mm));
    
    
    %% Perform the fit
    % We will search over the mRF transform parameters and lock the mRGC
    % transform parameters. Set upper and lower bounds on the mRF params
    % to be 1.25x the median values found across meridians (with a bit of
    % sign exponent trickery to handle the direction of negative params)
    lb = [rfInitialTransformParams./(1.25.^sign(rfInitialTransformParams)) rgcInitialTransformParams];
    ub = [rfInitialTransformParams.*(1.25.^sign(rfInitialTransformParams)) rgcInitialTransformParams];
    x0 = [rfInitialTransformParams rgcInitialTransformParams];
    
    % Set up the options
    options = optimoptions('fmincon', 'Display', 'none', 'ConstraintTolerance', 1);
    
    % Fit that sucker
    fitParams(mm,:) = fmincon(errorFunc,x0,[],[],[],[],lb,ub,nonlinconst,options);
    
    % Calculate and store displacement
    rgcDisplacementDegPolar(mm,:)=calcDisplacement(regularSupportPosDeg, mRGC_cumulative(fitParams(mm,:)), mRF_cumulative(fitParams(mm,:)));
    
    % Report the results for this meridian
    if p.Results.verbose
        zeroPoints = find(rgcDisplacementDegPolar(mm,:)==0);
        convergenceIdx = find(regularSupportPosDeg(zeroPoints) > 2,1);
        convergenceEccen = regularSupportPosDeg(zeroPoints(convergenceIdx));
        outLine = ['Polar angle: ' num2str(meridianAngles(mm)) ', max displacement: ' num2str(max(rgcDisplacementDegPolar(mm,:))) ', convergence eccen: ' num2str(convergenceEccen) '\n'];
        fprintf(outLine);
    end
    
    % Plot the results for this meridian
    if p.Results.makePlots
        % plot the displacement
        subplot(length(meridianAngles),2,mm*2);
        plot(regularSupportPosDeg,rgcDisplacementDegPolar(mm,:),'-r')
        axis off;
        ylim([-.5 3.0]);
        if mm == length(meridianAngles)
            axis on;
            xlabel('eccentricity [deg]');
            ylabel('RGC displacement [deg]');
        end
        
        % Plot the cumulative functions
        subplot(length(meridianAngles),2,mm*2-1);
        plot(regularSupportPosDeg,mRGC_cumulative(fitParams(mm,:)),'-k')
        axis off;
        if mm == length(meridianAngles)
            axis on;
            xlabel('eccentricity [deg]');
            ylabel('cells per sector');
        end
        hold on
        plot(regularSupportPosDeg,mRF_cumulative(fitParams(mm,:)),'-b')
        ylim([0 8e5]);
        hold off
        drawnow
    end
    
end % loop over meridians

% Create the displacement map
maxDisplacementDeg = max(rgcDisplacementDegPolar(:));
imP=rgcDisplacementDegPolar'./maxDisplacementDeg;
ImR = PolarToIm (imP, 0, 1, 1000, 1000);
displacementMapDeg = imrotate(ImR .* maxDisplacementDeg,-90);

% Plot the displacement map
if p.Results.makePlots
    figure
    climVals = [0,ceil(maxDisplacementDeg)];
    imagesc(displacementMapDeg, climVals);
    c = colorbar;
    c.Label.String='Radial RGC displacement [deg]';
    xlabel('Position [deg] temporal --> nasal');
    ylabel('Position [deg] inferior --> superior');
    k=linspace(-24,30,10);
    xticklabels(string(k))
    yticklabels(string(k))
end


end % calcDisplacementMap


%% LOCAL FUNCTIONS

function [c,ceq] = testRFGreaterThanRGC(regularSupportPosDeg, countPerRingRF, countPerRingRGC, displacementPointDeg)

% If there are any cumulative RGC values that are greater than the RF
% values at eccentricities less than the displacementPoint, then this
% violates the nonlinear fit constraint
withinRangeIdx = find(regularSupportPosDeg < displacementPointDeg);
c = sum(countPerRingRGC(withinRangeIdx) > countPerRingRF(withinRangeIdx));

ceq = []; % unused

end

function error = errorMatchingRFandRGC(regularSupportPosDeg, countPerRingRF, countPerRingRGC, displacementPointDeg)

% The error is calculated as the SSQ of the absolute difference between the
% cumulative RF and RGC counts at retinal positions past the point where
% displacement should have ended

withinRangeIdx = find(regularSupportPosDeg > displacementPointDeg);
error = sqrt(sum((countPerRingRGC(withinRangeIdx) - countPerRingRF(withinRangeIdx)).^2));
end

function displaceInDeg = calcDisplacement(regularSupportPosDeg, countPerRingRGC, countPerRingRF)

% Determine the sample resolution by a difference operation
tmp = diff(regularSupportPosDeg);
sampleResolutionDegrees = tmp(1);
% Measure the displacement (in degrees). First, for each cumulative RGC
% density value, identify the array index of the first
% cumulative RF density value that is equal or greater.
displaceInSamples=arrayfun(@(x) find(countPerRingRF>=x,1), countPerRingRGC,'UniformOutput',false);
% Now some array operations to get these values out of cells and in to a
% numeric vector
emptyCells = find(cellfun(@(x) isempty(x), displaceInSamples));
displaceInSamples(emptyCells)={NaN};
displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
% The displacement in array samples is the array index of each cumulative
% RGC density measurement, minus the array index of the matching RF density
% measurement. This difference is then multiplied by the sample resolution
% to get the displacement in degrees.
displaceInDeg = ((1:1:length(displaceInSamples))-displaceInSamples ) * sampleResolutionDegrees;
% Zero out negative values after the point of convergence
displaceInDeg(find(displaceInDeg < 0,1):end)=0;

end

