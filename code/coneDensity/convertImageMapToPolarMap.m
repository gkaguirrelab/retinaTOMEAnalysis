function polarMap = convertImageMapToPolarMap(imageMap, varargin)
% Returns a polar map with the same resolution as the
%
% Description:
%   Converts a square image to a polar map. This routine is a
%   wrapper for the function ImToPolar, created by Prakash Manandhar. The
%   wrapper handles a few properties of the conversion:
%     - The ImToPolar routine expects square matrix that has a
%       maximum value of unity.
%
% Inputs:
%   imageMap              - A square matrix of size imRdim that is the
%                           transform of the polarMap input.
%
% Optional key/value pairs:
%   imRdim                - The size of the imageMap to be returned in the
%                           X and Y dimension. If undefined, the imRdim
%                           will be set to preserve the full resolution of
%                           the input polarMap.
%
% Outputs:
%   polarMap              - An m x p matrix, consisting of m meridians and
%                           p eccentricity positions
%
% Examples:
%{
    % Test invertability of a sinusoidal grating
    
    % Create a grating (code taken from an Ione Fine class example)
    x=linspace(-pi, pi, 100);
    sf=6; % spatial freq in cycles per image
    sinewave=sin(x*sf);
    onematrix=ones(size(sinewave));
    imR=(onematrix'*sinewave);
    imRdim = 100;

    % Display the initial image
    figure
    subplot(1,4,1)
    imagesc(imR)
    colormap(gray)
    axis off;
    axis square;
    title('initial image');

    % Convert to polar and display
    imP = convertImageMapToPolarMap(imR);
    subplot(1,4,2)
    imagesc(imP)
    axis off;
    axis square;
    title('polar map');

    % Convert back to image and display
    imR2 = convertPolarMapToImageMap(imP,'imRdim',imRdim);
    subplot(1,4,3)
    imagesc(imR2)
    axis off;
    axis square;
    title('recovered image');

    % Show the error
    subplot(1,4,4)
    imagesc(imR-imR2)
    axis off;
    axis square;
    title('error');
%}    

%% Parse input and define variables
p = inputParser;

% Required
p.addRequired('polarMap',@isnumeric);

% Optional analysis params
p.addParameter('imRdim',size(imageMap,2)*2-1,@isnumeric);

% parse
p.parse(imageMap,varargin{:})

% Pull out the result from the param parser
imRdim = p.Results.imRdim;

% Rotate the imageMap to correspond to what is expected for ImToPolar
imageMap = imrotate(imageMap ,90);

% Scale the imageMap to have a max value of unity, saving the scaler
maxMapDensity = max(imageMap(:));
imR=imageMap./maxMapDensity;

% Convert from polar to image
imP = ImToPolar (imR, 0, 1, imRdim, imRdim);

% Re-scale to original max
polarMap = imP' .* maxMapDensity;

end