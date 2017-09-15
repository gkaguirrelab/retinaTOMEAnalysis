function [ mRGCDensitySqDeg, midgetFraction ] = transformRGCToMidgetRGCDensity( regularSupportPosDeg, rgcDensitySqDeg, varargin )
% transformRGCToMidgetRGCDensity( regularSupportPosDeg, rgcDensitySqDeg, varargin )
%
% This function transforms a vector of RGC density measurements into
% corresponding midget RGC density measurements. The eccentricity location
% of the measurements are not used explicitly in the calculation, instead
% the cumulative sum of the RGC density.
%
% This transformation is under the control of four parameters. The first is
% the f0 value (from Eq 7 of Watson 2014 JoV) that defines the midget
% fraction at the fovea. The next three are parameters of a reciprocal
% function.
%
% If x is the proportion of the cumulative RGC density function at a given
% location, then:
%
%   midgetFraction =  f0 - [(1./(a+(b.* log10(x) )))+c]
%
% where f0 is the maximum midget fraction found at the fovea, given by
% Watson / Drasdo as 0.8928.
%
% INPUT
%   regularSupportPosDeg - The eccentricity values that support
%       rgcDensitySqDeg.
%   rgcDensitySqDeg - RGC density at each eccentricity location in units of
%      cells / degree^2
%
% OUTPUT
%   mRGCDensitySqDeg - midget RGC density at each eccentricity location
%   midgetFraction - the midget fraction at each eccentricity location
%
% OPTIONS
%   referenceEccen - the reference eccentricity for the proportion of
%       the cumulative RGC density. The proportion function will have a
%       value of unity at this point.
%   watsonEq8_f0 - The midget fraction assigned to the fovea
%   recipFitParams - parameters of the reciprocal fit the defines the
%      transformation
%   verbose - Controls text output to console
%   makePlots - Do we make a figure?

%% Parse input and define variables
p = inputParser;

% required input
p.addRequired('regularSupportPosDeg',@isnumeric);
p.addRequired('rgcDensitySqDeg',@isnumeric);

% Optional anaysis params
p.addParameter('referenceEccen',30,@isnumeric);
p.addParameter('watsonEq8_f0',0.8928,@isnumeric);
p.addParameter('recipFitParams',[2.4026 -8.0877 -0.0139],@isnumeric);

% parse
p.parse(regularSupportPosDeg, rgcDensitySqDeg, varargin{:})


%% House keeping and setup

% Define a three-parameter reciprocal function that will be used to fit the
% modeled relationship
recipFunc = fittype('(1./(a+(b.*x)))+c','independent','x','dependent','y');

% Obtain the cumulative RGC function
RGC_ringcount = calcCumulative(regularSupportPosDeg,rgcDensitySqDeg);

% Find the index position in the regularSupportPosDeg that is as close
% as possible to the referenceEccen
[ ~, refPointIdx ] = min(abs(regularSupportPosDeg-p.Results.referenceEccen));
% Calculate a proportion of the cumulative RGC density counts, relative
% to the reference point (which is assigned a value of unity)
propRGC_ringcount=RGC_ringcount./RGC_ringcount(refPointIdx);

% Because we are going to be working with a log transform, set any zero
% proportion values to the minimum, non-zero proportion value
zeroPoints=find(propRGC_ringcount==0);
if ~isempty(zeroPoints)
    propRGC_ringcount(zeroPoints)=min(propRGC_ringcount(find(propRGC_ringcount~=0)));
end

% Calculate the midgetFraction based upon the propRGC_ringcount
midgetFraction = p.Results.watsonEq8_f0-recipFunc(p.Results.recipFitParams(1),p.Results.recipFitParams(2),p.Results.recipFitParams(3),log10(propRGC_ringcount));

% Scale the rgcDensity by the midget fraction
mRGCDensitySqDeg = rgcDensitySqDeg .* midgetFraction;

% NaNs can happen; set them to zero
badIdx = find(isnan(mRGCDensitySqDeg));
if ~isempty(badIdx)
    mRGCDensitySqDeg(badIdx)=0;
end

end % function


