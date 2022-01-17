function [smoothMap] = fillAndSmoothMap(Map,sampleBaseX,sampleBaseY,varargin)
% Interpolate and smooth the values in a map
%
% Description:
%   Needed.
%
% Inputs:
%   Map                   - 2D matrix that the interp and smoothing are
%                           applied to.
%   sampleBaseX           - 2D support fot the X dim.
%   sampleBaseY           - 2D support fot the Y dim.
%
% Optional key/value pairs:
%  'interpMethod'         - Type of interpolation to remove holes in mesh.
%                           Legal values options:
%                              'linear'	- Triangulation-based linear 
%                                       interpolation (default) supporting
%                                       2-D and 3-D interpolation.
%                               'nearest' - Triangulation-based nearest 
%                                       neighbor interpolation supporting
%                                       2-D and 3-D interpolation.
%                                       Discontinuous
%                               'natural' - Triangulation-based natural 
%                                       neighbor interpolation supporting
%                                       2-D and 3-D interpolation. This
%                                       method is an efficient tradeoff
%                                       between linear and cubic.
%                               'cubic' - Triangulation-based cubic 
%                                       interpolation supporting 2-D
%                                       interpolation only.
%                               'v4'  - Biharmonic spline interpolation. 
%                                       Unlike the other methods, this
%                                       interpolation is not based on a
%                                       triangulation.
%  'kernelType'           - Kernel to be used for smoothing. Options:
%                               'gaussian' 
%                               'disk'
%                               'average'
%  'sigma'                - Standard deviation of smoothing kernel
%  'hsize'                - Size of the smoothing kernel.
%  'verbose'              - say stuff
%
% Outputs:
%   smoothMap             - The smoothed map
%
% Example call:
%   function [smoothMap] = fillAndSmoothMap(Map,sampleBaseX,sampleBaseY,varargin)
%
% mab 2017
%


p = inputParser;

% Optional anaysis params
p.addParameter('interpMethod','linear');
p.addParameter('kernelType','gaussian');
p.addParameter('sigma',3,@isnumeric);
p.addParameter('hsize',15,@isnumeric);


% Optional display params
p.addParameter('verbose',false,@islogical);

% parse
p.parse(varargin{:})


%% fix tears in warped image
%// identify indices valid for the map
idxgood=~(isnan(sampleBaseX) | isnan(sampleBaseY) | isnan(Map)); 

%// re-interpolate scattered data over the sampleBase grid
map_filled = griddata( sampleBaseX(idxgood),sampleBaseY(idxgood),Map(idxgood), sampleBaseX, sampleBaseY); 


%% Smooth the final map
switch p.Results.kernelType
    case 'gaussian'
        H = fspecial(p.Results.kernelType,p.Results.hsize,p.Results.sigma);
    case 'disk'
        H = fspecial(p.Results.kernelType,round(p.Results.hsize)/2);
    case 'average'
        H = fspecial(p.Results.kernelType,p.Results.hsize);
end
       
smoothMap = imfilter(map_filled,H,'replicate'); 

if p.Results.verbose
    figure;
    imgacesc(smoothMap)
end

end