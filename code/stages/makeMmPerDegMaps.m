function makeMmPerDegMaps(saveDir, varargin)
% Creates and saves maps of the conversion of degrees to mm on the retina
%
% Description:
%   Each position in the visual field projects to a position on the retina.
%   Position in the visual field is expressed in degrees of visual angle.
%   This routine calculates, for any position in the visual field, the mm
%   of retina that are subtended by a degree of visual angle at that
%   position. A general rule of thumb is that one mm of retina corresponds
%   to three degrees of visual angle. The exact value will vary principally
%   based upon the axial length of the eye, and to a lesser extent upon the
%   position in the visual field as influenced by the optics of the cornea,
%   the misalignment of the optical and visual axes of the eye, and the
%   particulars of the radii of the ellipsoidal vitreous chamber.
%
%   This routine loads the eye biometric parameters measured for the TOME
%   subjects and uses these to create a model eye for each subject. Ray
%   tracing in that model eye is then conducted to produce a map of the mm
%   of retina per degree of visual angle on the retinal surface. This map
%   is then fit with a polynomial surface and the parameters retained. To
%   obtain the mmPerDeg for a given subject at a given position in the
%   visual field, relative to the fovea, use:
%
%       mmPerDeg = polyval(mmPerDegPolyFit{ss},[horizDeg, vertDeg])
%
%   The position on the retina will be the negative of these.  
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('saveDir',@ischar);

% Optional analysis params
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);

%% Parse and check the parameters
p.parse(saveDir, varargin{:});

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% Create an empty cell array to hold the results
mmPerDegPolyFit = {};

parfor ii = 1:size(subjectTable.AOSO_ID,1)
%for ii = 1:length(subjectTable.AOSO_ID)

    % Assemble the corneal curvature, spherical error, and axial length.
    axialLength = subjectTable.Axial_Length_average(ii);
    SR = subjectTable.Spherical_Error_average(ii);
    cc = subjectTable.modelEyeK_OD(ii);
    cc=str2num(cc{1});

    % Create the model eye
    eye = modelEyeParameters('axialLength',axialLength,'sphericalAmetropia',SR,'measuredCornealCurvature',cc,'calcLandmarkFovea',true);
    
    % Define the visual field domain over which we will make the measure
    horizVals = -30:15:30;
    vertVals = -30:15:30;

    % Define an empty matrix to hold the results
    mmPerDeg = nan(length(horizVals),length(vertVals));

    % Define the delta deg
    deltaDegEuclidean = 1;
    deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];
    
    % Loop over horizontal and vertical field positions
    for jj = 1:length(horizVals)
        for kk = 1:length(vertVals)
            % The position in the field relative to the optical axis of the
            % eye
            degField = [horizVals(jj) vertVals(kk)] + eye.landmarks.fovea.degField(1:2);
            % Obtain the retinal points that are delta degrees on either
            % side of the specified degree field position
            [~,X0] = calcRetinaFieldPoint( eye, degField - deltaAngles./2);
            [~,X1] = calcRetinaFieldPoint( eye, degField + deltaAngles./2);
            mmPerDeg(jj,kk) = norm(X0-X1) / norm(deltaAngles);
        end
    end
    % Fit a polynomial surface to the measure
    [X,Y]=meshgrid(horizVals,vertVals);
    pp = fit([X(:),Y(:)],mmPerDeg(:),'poly33');
    % Give some console update
    fprintf(['Done subject ' num2str(subjectTable.AOSO_ID(ii)) '\n']);
    % Store the map in a cell array that accumulates across the parfor
    mmPerDegPolyFit(ii) = {pp};    
end

% Write out the results
outfile = fullfile(saveDir,'mmPerDegPolyFit.mat');
save(outfile,'mmPerDegPolyFit');

end % Main


