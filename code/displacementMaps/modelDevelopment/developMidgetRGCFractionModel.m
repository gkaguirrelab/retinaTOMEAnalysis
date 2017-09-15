function [ fitParams, figHandle ] = developMidgetRGCFractionModel( varargin )
% developMidgetRGCFractionModel( varargin )
%
% Our goal is to be able to derive the fraction of RGCs that are midget
% RGCs at a point in the retina, but without reference to the polar angle
% or eccentricity of the measurement.
%
% To solve this problem, we first take the cumulative of the RGC density
% function along a meridian. We then express this function as a proportion
% of the cumuluative density value at any point to the cumulative denisty
% value at a specific eccentricity location (e.g., at 30 degrees). Thus,
% the proportional cumulative density of RGCs is equal to unity (1) at the
% reference eccentricity. The value of the reference eccentricity should be
% sufficiently large so that the RGC density function is in the
% monotonically decreasing phase.
%
% Next, the log10 of the proportional, cumulative RGC density function is
% related to the predicted midget fraction for that location, as derived by
% eq 7 of Watson 2014 JOV (which itself is taken from Drasdo 2007). We
% observe that a three-component reciprocal function fits this relationship
% well.
%
% We obtain the fit parameters for this relationship for each of the four
% meridians for which we have data, and take the median fit parameters
% across merdians. This gives us a functional form that can be used to
% convert RGC density to midget RGC density knowing only the proportion of
% the cumulative RGC density function at the location for which the
% measurement was made. If we term this proportion value x, then the
% relationship is:
%
%   midgetFraction =  f0 - [(1./(a+(b.* log10(x) )))+c]
%
% where f0 is the maximum midget fraction found at the fovea, given by
% Watson / Drasdo as 0.8928.
%
% OUTPUT
%   fitParams - The median values of [a b c] across the four merdians
%   figHandle - Handle to a figure showing the fits. Empty if no plotting
%      requested.
%
% OPTIONS
%   referenceEccen - the reference eccentricity for the proportion of
%       the cumulative RGC density. The proportion function will have a
%       value of unity at this point.
%   supportResolutionDegrees - the resolution (in degrees) for which the
%       calculations are made.
%   supportEccenMaxDegrees - the maximum eccentricity used for modeling
%   meridianNames - Cell array of the text string names of the meridia
%   meridianAngles - Polar angle values assigned to the meridians
%   watsonEq8_f0 - The midget fraction assigned to the fovea
%   watsonEq8_rm - decay parameter controlling the midget fraction function
%   recipFitStartPoint - initial values used for the fit of the
%       three-parameter reciprocal function. It is critical that the second
%       parameter (b) has a negative start point to fit this inverted
%       reciprocal.
%   verbose - Controls text output to console
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% required input
%p.addRequired('grayVideoName',@isstr);

% Optional anaysis params
p.addParameter('referenceEccen',30,@isnumeric);
p.addParameter('supportResolutionDegrees',0.01,@isnumeric);
p.addParameter('supportEccenMaxDegrees',30,@isnumeric);
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@cell);
p.addParameter('watsonEq8_f0',0.8928,@isnumeric);
p.addParameter('watsonEq8_rm',41.03,@isnumeric);
p.addParameter('recipFitStartPoint',[3 -8 0],@isnumeric);

% Optional display params
p.addParameter('verbose',true,@islogical);
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(varargin{:})


%% House keeping and setup

% Define a three-parameter reciprocal function that will be used to fit the
% modeled relationship
recipFunc = fittype('(1./(a+(b.*x)))+c','independent','x','dependent','y');

% Define the regular-spaced eccentricity support over which we will model
% the anatomical retinal functions
regularSupportPosDeg = ...
    0:p.Results.supportResolutionDegrees:p.Results.supportEccenMaxDegrees;

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
else
    figHandle=[];
end

%% Loop over the meridians

for mm = 1:length(p.Results.meridianAngles)
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ RGCDensitySqDeg, nativeSupportPosDeg ] = getCurcioRGCDensityByEccen( p.Results.meridianAngles(mm) );
    
    % remove nan values
    isvalididx=find(~isnan(RGCDensitySqDeg)  );
    nativeSupportPosDeg = nativeSupportPosDeg(isvalididx);
    RGCDensitySqDeg = RGCDensitySqDeg(isvalididx);
    
    % Fit a spline to the RGC density data
    RGCDensityFit = fit(nativeSupportPosDeg,RGCDensitySqDeg,'smoothingspline', 'Exclude',find(isnan(RGCDensitySqDeg)),'SmoothingParam', 1);
    
    % Obtain the cumulative RGC function
    RGC_ringcount = calcCumulative(regularSupportPosDeg,RGCDensityFit(regularSupportPosDeg)');
    
    % Find the index position in the regularSupportPosDeg that is as close
    % as possible to the referenceEccen
    refPointIdx= ...
        find((regularSupportPosDeg-p.Results.referenceEccen)== ...
        min(abs(regularSupportPosDeg-p.Results.referenceEccen)));
    % Calculate a proportion of the cumulative RGC density counts, relative
    % to the reference point (which is assigned a value of unity)
    propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);
    
    % Because we are going to be working with a log transform, set any zero
    % proportion values to the minimum, non-zero proportion value
    zeroPoints=find(propRGC_ringcount==0);
    if ~isempty(zeroPoints)
        propRGC_ringcount(zeroPoints)=min(propRGC_ringcount(find(propRGC_ringcount~=0)));
    end
    
    % Obtain the Watson midget fraction as a function of eccentricity
    midgetFractionByEccen = calcWatsonMidgetFractionByEccen(regularSupportPosDeg,p.Results.watsonEq8_f0,p.Results.watsonEq8_rm);
    
    % Fit the reciprocal model that relates log10(proportionRGC) to midget
    % fraction. First, define a weight function to lock the f0 value
    weights=ones(1,length(propRGC_ringcount));
    weights(1)=1e5;
    
    % Perform the fit and save the param values.
    % The x values are (-1) * log10(proportionRGC)
    % The y values are the f0 value minus the midget fraction
    
    % turn off some warnings that are produced by the fit and plot
    warningState = warning;
    warning('off','curvefit:fit:complexXusingOnlyReal');
    warning('off','MATLAB:plot:IgnoreImaginaryXYPart');

    recipFit= ...
        fit(log10(propRGC_ringcount'),p.Results.watsonEq8_f0-midgetFractionByEccen', ...
        recipFunc,'Weights',weights,'StartPoint',p.Results.recipFitStartPoint );
    fitParams(mm,1)=recipFit.a;
    fitParams(mm,2)=recipFit.b;
    fitParams(mm,3)=recipFit.c;
    
    % Add the data to the figure
    if p.Results.makePlots
        plot(log10(propRGC_ringcount),midgetFractionByEccen,p.Results.meridianSymbols{mm},'color',[.8 .8 .8])
        hold on
        ylim([0.4 1]);
        xlim([-12 1]);
        xlabel('log10 proportion cumulative RGC density count');
        ylabel('midget fraction');
    end % if we are plotting
    
    % restore the saved warning state
    warning(warningState);
    
end % loop over meridians

% Now calculate the median parm values across the meridians and add a fit
% line to the plot
fitParams=median(fitParams);

if p.Results.makePlots
    xFit=logspace(-12,0,100);
    plot( log10(xFit), ...
        p.Results.watsonEq8_f0-recipFunc(fitParams(1),fitParams(2),fitParams(3),log10(xFit)),'-r')
    legend({p.Results.meridianNames{:} 'fit'},'Location','southwest');
    title('midget fraction as a function of relative RGC cumulative density');
    drawnow
end


end % function


