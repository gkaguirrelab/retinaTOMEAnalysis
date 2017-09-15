function [ fitParams, figHandle ] = developMidgetRFFractionModel( varargin )
% developMidgetRFFractionModel( varargin )
%
% Examine the relationship between cone density and midget
% receptive field density. Note that the mRF density reflects both the
% ratio of convergence of cones onto retinal ganglion cells, as well as the
% fraction of retinal ganglion cells at any one point that are midgets.
%
% Judging from the work of Watson (2014), the ratio at the fovea, would be
% expected to be 1.786 : 1, as there is believed to be a 2:1 cone:RGC
% ratio, and the midget fraction is modeled as being 0.8928. However, a
% higher value (1.9) was observed to best fit the data across meridians.
%
% For each meridian, We first load the Curcio cone density measurements.
% At each of the sampled eccentricity values, we obtain the midgetRF
% density as provided by Eq 8 of Watson 2014.
%
% We observe that the ratio of midgetRF density to cone density as a
% function of log10 cone density is sigmoidal. To support generalization
% across datasets, we express the x-axis as log10 of the proportion of
% maximum cone density observed in the fovea. We then fit this relationship
% with the function:
%
%	mRFtoConeDensityRatio = minRatio+(maxRatio-minRatio)./(1+(x./inflect).^slope)
%
% where maxRatio and minRatio are locked at the values of 1.9 and 0, respectively.
%
% We calculate the fit of this function across each of the four meridians,
% then take the median of the two free fit parameters (slope and inflect)
% and return these.
%
%
% OUTPUT
%   fitParams - The median values of [slope inflect] across the four merdians
%   figHandle - Handle to a figure showing the fits. Empty if no plotting
%      requested.
%
% OPTIONS
%   supportEccenMaxDegrees - the maximum eccentricity used for modeling.
%       This value should be sufficiently high so that we are in the
%       asymptote range of cone density
%   meridianNames - Cell array of the text string names of the meridia
%   meridianAngles - Polar angle values assigned to the meridians
%   maxConeDensity - The maximum cone density at the fovea (couts / deg^2).
%       The default value is from Curcio 1990. If set to empty, the maximum
%       value from coneDensitySqDeg is used.
%   minRatio - The minimum value of the mRF:cone density ratio. Set to zero
%       as the functions appear to asymptote close to this value.
%   maxRatio - The maximuim value of the mRF:cone density ratio. Set to
%       1.9 following the logic outlined above.
%   logitFitStartPoint - initial values used for the slope and inflection
%       point parameters of the logisic fit. Informed by examination of
%       typical data.
%   verbose - Controls text output to console
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% required input
%p.addRequired('grayVideoName',@isstr);

% Optional anaysis params
p.addParameter('referenceEccen',30,@isnumeric);
p.addParameter('supportEccenMaxDegrees',60,@isnumeric);
p.addParameter('meridianNames',{'Nasal' 'Superior' 'Temporal' 'Inferior'},@iscell);
p.addParameter('meridianAngles',[0, 90, 180, 270],@isnumeric);
p.addParameter('meridianSymbols',{'.','x','o','^'},@cell);
p.addParameter('maxConeDensity',1.4806e+04,@(x)(isempty(x) | isnumeric(x)));
p.addParameter('minRatio',0,@isnumeric);
p.addParameter('maxRatio',1.9,@isnumeric);
p.addParameter('logitFitStartPoint',[3,-1],@isnumeric);

% Optional display params
p.addParameter('verbose',true,@islogical);
p.addParameter('makePlots',true,@islogical);

% parse
p.parse(varargin{:})

% Set the maxConeDensity value
if isempty(p.Results.maxConeDensity)
    maxConeDensity = max(coneDensitySqDeg);
else
    maxConeDensity = p.Results.maxConeDensity;
end

%% House keeping and setup

% Define a four-parameter logistic function that will be used to fit the
% modeled relationship. Two of the parameters (max and min asymptote) are
% locked by the passed parameter
logisticFunc = fittype( @(slope,inflect,minRatio,maxRatio,x) minRatio+(maxRatio-minRatio)./(1+sign(x./inflect).*abs((x./inflect).^slope)), ...
    'independent','x','dependent','y','problem',{'minRatio','maxRatio'});

% Prepare a figure if requested
if p.Results.makePlots
    figHandle=figure;
else
    figHandle=[];
end

%% Loop over the meridians

for mm = 1:length(p.Results.meridianAngles)
    
    % Load the Cone density Data from Curcio 1990:
    [nativeSupportPosDeg,coneDensitySqDeg] = getCurcioConeDensitySqDeg(p.Results.meridianAngles(mm));
    
    % remove nan values
    isvalididx=find(~isnan(coneDensitySqDeg)  );
    nativeSupportPosDeg = nativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    
    % calculate the mRF density at these eccentricity locations using
    % Watson equation 8.
    [ midgetRFDensitySqDeg ] = getWatsonMidgetReceptiveFieldDensityByEccen(nativeSupportPosDeg, p.Results.meridianAngles(mm));
    
    % Remove nans and points beyond the modeled eccentricity bound
    isvalididx=find(~isnan(midgetRFDensitySqDeg).*~isnan(coneDensitySqDeg) .* (nativeSupportPosDeg < p.Results.supportEccenMaxDegrees) );
    nativeSupportPosDeg = nativeSupportPosDeg(isvalididx);
    coneDensitySqDeg = coneDensitySqDeg(isvalididx);
    midgetRFDensitySqDeg = midgetRFDensitySqDeg(isvalididx);
    
    % Define the ratio function.
    midgetRFtoConeRatio = (midgetRFDensitySqDeg ./ coneDensitySqDeg)';
    
    % Define the x-axis as the log10 of the proportion of max cone density
    x = log10(coneDensitySqDeg ./ maxConeDensity)';
        
    % Perform the logistic fit. Note that the max and min asymptote are
    % pinned by the passed parameters
    logitFit = fit(x,midgetRFtoConeRatio,logisticFunc, ...
        'problem',{p.Results.minRatio, p.Results.maxRatio}, ...
        'StartPoint',p.Results.logitFitStartPoint, ...
        'Lower',[0,-2],'Upper',[10,0] );
    
    % Distribute the logitFit params to a loop variable
    fitParams(mm,1)=logitFit.slope;
    fitParams(mm,2)=logitFit.inflect;
    
    % Add the data to the figure
    if p.Results.makePlots
        plot(x,midgetRFtoConeRatio,p.Results.meridianSymbols{mm},'Color',[.7 .7 .7]);
        hold on
        xlabel('log10 proportion max cone density');
        ylabel('mRF density : cone density');
    end % if we are plotting
    
end % loop over meridians

% Now calculate the median parm values across the meridians and add a fit
% line to the plot
fitParams=median(fitParams);

if p.Results.makePlots
    xFit= -2:.01:0;
    plot( xFit, ...
        logisticFunc(fitParams(1), fitParams(2), p.Results.minRatio, p.Results.maxRatio, xFit),'-r')
    legend({p.Results.meridianNames{:} 'fit'},'Location','southeast');
    title('midget RF : cone ratio as a function of relative cone density');
    drawnow
end


end % function


