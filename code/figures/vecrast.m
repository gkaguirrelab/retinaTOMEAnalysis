function vecrast(figureHandle, filename, resolution, stack, exportType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodoros Michelis, 6 October 2017
% TUDelft, Aerospace Engineering, Aerodynamics
% t.michelis@tudelft.nl
%
%
% D E S C R I P T I O N:
% vecrast is a function that allows to automatically save a figure with
% mixed vector and raster content. More specifically, two copies of the
% figure of interest are created, rasterFigure and vectorFigure. Patches,
% surfaces, contours, images, and lights are kept in rasterFigure but
% removed from vectorFigure. rasterFigure is then saved as a temporary
% .png image with the required resolution. The .png file is subsequently
% inserted into the vectorFigure, and the result is saved in a single
% vector file.
%
%
% I N P U T:
% vecrast(figureHandle, filename, resolution, stack, exportType)
%   figureHandle:   Handle of the figure of interest
%   filename:       Baseline name string of output file WITHOUT the extension.
%   resolution:     Desired resolution of rasterising in dpi
%   stack:          'top' or 'bottom'. Stacking of raster image with
%                       respect to axis in vector figure, see examples below.
%   exportType:     'pdf' or 'eps'. Export file type for the output file.
%
%
% N O T E S:
% - The graphics smoothing (anti-aliasing) is turned off for the raster
%   figure. This improves sharpness at the borders of the image and at the
%   same time greatly reduces file size. You may change this option in the
%   script by setting 'GraphicsSmoothing', 'on' (line 84).
% - A resolution of no less than 300 dpi is advised. This ensures that
%   interpolation at the edges of the raster image does not cause the image
%   to bleed outside the prescribed axis (make a test with 20dpi on the
%   first example and you will see what I mean).
% - The stacking option has been introduced to accomodate 2D and 3D plots
%   which require the image behind or in front the axis, respectively. This
%   difference can be seen in the examples below.
% - I strongly advise to take a look at the tightPlots function that allows
%   setting exact sizes of figures.

% E X A M P L E   1:
%   clear all; close all; clc;
%   Z = peaks(20);
%   contourf(Z,10)
%   vecrast(gcf, 'example1', 300, 'bottom', 'pdf');

% E X A M P L E   2:
%   clear all; close all; clc;
%   [X,Y] = meshgrid(1:0.4:10, 1:0.4:20);
%   Z = sin(X) + cos(Y);
%   surf(X,Y,Z)
%   vecrast(gcf, 'example2', 300, 'top', 'pdf');

% Thanks to Jonathan Kohler, Kerry and Bob DA for their input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some checks of the input ------------------------------------------------
if strcmp(stack, 'top') + strcmp(stack, 'bottom') == 0
    error('Stack must be ''top'' or ''bottom''');
end
if strcmp(exportType, 'pdf') + strcmp(exportType, 'eps') == 0
    error('Stack must be ''pdf'' or ''eps''');
end

% Ensure figure has finished drawing
drawnow;

% Set figure units to points
set(figureHandle, 'units', 'points');

% Ensure figure size and paper size are the same
figurePosition = get(figureHandle, 'Position');
set(figureHandle, 'PaperUnits', 'points', 'PaperSize', [figurePosition(3) figurePosition(4)])
set(figureHandle, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 figurePosition(3) figurePosition(4)]);

% Create a copy of the figure and remove smoothness in raster figure
rasterFigure = copyobj(figureHandle, groot);
vectorFigure = copyobj(figureHandle, groot);
set(rasterFigure, 'GraphicsSmoothing', 'off', 'color', 'w');
set(vectorFigure, 'GraphicsSmoothing', 'off', 'color', 'w');

% Fix vector image axis limits based on the original figure
% (this step is necessary if these limits have not been defined)
axesHandle = findall(vectorFigure, 'type', 'axes');
for ii = 1:length(axesHandle)
    xlim(axesHandle(ii), 'manual');
    ylim(axesHandle(ii), 'manual');
    zlim(axesHandle(ii), 'manual');
end

% Create axis in vector figure to fill with raster image
rasterAxis = axes(vectorFigure, 'color', 'none', 'box', 'off', 'units', 'points');
set(rasterAxis, 'position', [0 0 figurePosition(3) figurePosition(4)]);
uistack(rasterAxis, stack);

% Ensure fontsizes are the same in all figures
figText = findall(figureHandle, 'type', 'text');
rastText = findall(rasterFigure, 'type', 'text');
vecText = findall(vectorFigure, 'type', 'text');
for ii=1:length(figText)
    set(rastText(ii), 'FontSize', get(figText(ii), 'FontSize'));
    set(vecText(ii), 'FontSize', get(figText(ii), 'FontSize'));
end

% Raster Figure ----------------------------------------------------------
% Select what to remove from raster figure
axesHandle = findall(rasterFigure, 'type', 'axes');
axesPosition = get(axesHandle,'position'); % Fix: get axes size
set(axesHandle, 'color', 'none');
for ii = 1:length(axesHandle)
    contents = findall(axesHandle(ii));
    toKeep = [...
        findall(axesHandle(ii), 'type', 'patch');...
        findall(axesHandle(ii), 'type', 'surface');...
        findall(axesHandle(ii), 'type', 'contour')...
        findall(axesHandle(ii), 'type', 'image');...
        findall(axesHandle(ii), 'type', 'light')];
    toRemove = setxor(contents, toKeep);
    set(toRemove, 'visible', 'off');
end
if ~isa(axesHandle,'matlab.graphics.axis.Axes')
    for ii = 1:length(axesHandle)
        set(axesHandle(ii),'position',axesPosition{ii}); % Fix: restore original axes size
    end
end

% Remove all annotations from raster figure
annotations = findall(rasterFigure, 'Tag', 'scribeOverlay');
for ii = 1:length(annotations)
    set(annotations(ii), 'visible', 'off');
end

% Hide all colorbars and legends from raster figure
colorbarHandle = findall(rasterFigure, 'type', 'colorbar');
legendHandle = findall(rasterFigure, 'type', 'legend');
set([colorbarHandle; legendHandle], 'visible', 'off');

% Print rasterFigure on temporary .png
% ('-loose' ensures that the bounding box of the figure is not cropped)
print(rasterFigure, [filename 'Temp.png'], '-dpng', ['-r' num2str(resolution) ], '-loose', '-opengl');
close(rasterFigure);

% Vector Figure -----------------------------------------------------------
% Select what to keep in vector figure
axesHandle = findall(vectorFigure, 'type', 'axes');
set(axesHandle, 'color', 'none');
for ii = 1:length(axesHandle)
    toRemove = [...
        findall(axesHandle(ii), 'type', 'patch');...
        findall(axesHandle(ii), 'type', 'surface');...
        findall(axesHandle(ii), 'type', 'contour');...
        findall(axesHandle(ii), 'type', 'image');...
        findall(axesHandle(ii), 'type', 'light');...
        ];
    set(toRemove, 'visible', 'off');
end

% Insert Raster image into the vector figure
[A, ~, alpha] = imread([filename 'Temp.png']);

if isempty(alpha)==1
    imagesc(rasterAxis, A);
else
    imagesc(rasterAxis, A, 'alphaData', alpha);
end
axis(rasterAxis, 'off');

% Ensure that all colorbar ticks match the original figure
hcbOriginal = findall(figureHandle, 'type', 'colorbar');
hcbVector = findall(vectorFigure, 'type', 'colorbar');
if ~isempty(hcbOriginal)
    cbLimits = hcbOriginal.Limits;
    hcbVector.Ticks = (hcbOriginal.Ticks - cbLimits(1))/diff(cbLimits);
    hcbVector.TickLabels = hcbOriginal.TickLabels;
end

% Bring all annotations on top
annotations = findall(vectorFigure, 'Tag', 'scribeOverlay');
for ii = 1:length(annotations)
    uistack(annotations(ii), 'top');
end
% Ensure figure has finished drawing
drawnow;

% Finalise ----------------------------------------------------------------
% Remove raster image from directory
delete([filename 'Temp.png']); % COMMENT THIS IF YOU WANT TO KEEP PNG FILE

% Print and close the combined vector-raster figure
% set(vectorFigure, 'Renderer', 'painters');
if strcmp(exportType, 'pdf') == 1
    print(vectorFigure, [filename '.pdf'], '-dpdf', '-loose', '-painters', ['-r' num2str(resolution) ]);
elseif strcmp(exportType, 'eps') == 1
    print(vectorFigure, [filename '.eps'], '-depsc2', '-loose', '-painters');
end

close(vectorFigure);

end