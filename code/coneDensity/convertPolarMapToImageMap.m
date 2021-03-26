function imageMap = convertPolarMapToImageMap(polarMap, varargin)
% Returns an image map with the same resolution as the polar map
%
% Description:
%   Converts the polarMap matrix to a square image. This routine is a
%   wrapper for the function PolarToIm, created by Prakash Manandhar. The
%   wrapper handles a few properties of the conversion:
%     - The PolarToIm routine expects an input polar matrix that has a
%       maximum value of unity.
%     - The resulting image must be rotated to match our image convention
%     - We set points beyond the circular defined boudary to nan
%
% Inputs:
%   polarMap              - An m x p matrix, consisting of m meridians and
%                           p eccentricity positions
% Optional key/value pairs:
%   imRdim                - The size of the imageMap to be returned in the
%                           X and Y dimension. If undefined, the imRdim
%                           will be set to preserve the full resolution of
%                           the input polarMap.
%
% Outputs:
%   imageMap              - A square matrix of size imRdim that is the
%                           transform of the polarMap input.
%

%% Parse input and define variables
p = inputParser;

% Required
p.addRequired('polarMap',@isnumeric);

% Optional analysis params
p.addParameter('imRdim',size(polarMap,2)*2-1,@isnumeric);

% parse
p.parse(polarMap,varargin{:})

% Pull out the result from the param parser
imRdim = p.Results.imRdim;

% Scale the polar map to have a max value of unity, saving the scaler
maxMapDensity = max(polarMap(:));
imP=polarMap'./maxMapDensity;

% Convert from polar to image
imR = PolarToIm (imP, 0, 1, imRdim, imRdim);

% Fix the rotation of the map and re-scale to original max
imageMap = imrotate(imR .* maxMapDensity,-90);

% Set values beyond the radius of the provided meridians to nans
distanceMap = zeros(size(imageMap));
distanceMap(ceil(imRdim/2),ceil(imRdim/2))=1;
distanceMap=bwdist(distanceMap,'euclidean');
imageMap(distanceMap>(floor(imRdim/2)))=nan;

end