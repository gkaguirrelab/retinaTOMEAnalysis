function makeMmPerDegMaps(saveDir, varargin)
% Creates and saves maps of the conversion of degrees to mm on the retina
%
% Description:
%   Each position in the visual field projects to a position on the retina.
%   Position in the visual field is expressed in degrees of visual angle.
%   This routine calculates, for any position in the visual field, the mm
%   of retina that are subtended by a degree of visual angle at that
%   position. A general rule of thumb is that one degree of visual angle
%   equals roughly three millimeters. The exact value will vary principally
%   based upon the axial length of the eye, and to a lesser extent upon the
%   position in the visual field as influenced by the optics of the cornea,
%   the misalignment of the optical and visual axes of the eye, and the
%   particulars of the radii of the ellipsoidal vitreous chamber.
%
%   This routine loads the eye biometric parameters measured for the TOME
%   subjects and uses these to create a model eye for each subject. Ray
%   tracing in that model eye is then conducted to produce a map of the mm
%   of retina per degree of visual angle on the retinal surface.
%
%   The resulting map is 13 x 13 samples, and has the following properties:
%       mmPerDeg(1,1) -- nasal, superior
%       mmPerDeg(13,1) -- temporal, superior
%       mmPerDeg(1,13) -- nasal, inferior
%       mmPerDeg(13,13) -- temporal, inferior
%       mmPerDeg(7,7) -- fovea
%
%   The model assumes that all subjects have the same alpha angles (the
%   angle between the visual and optical axes of the eye). While this is an
%   inaccurate assumption, we lack this measurement for each subject, and
%   instead use the mean value for all subjects.
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('saveDir',@ischar);

% Optional analysis params
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);
p.addParameter('alpha',[5.45 2.5 0],@isnumeric);

%% Parse and check the parameters
p.parse(saveDir, varargin{:});


% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);

% Create an empty cell array to hold the results
resultSet = {};

parfor (ii = 1:length(subjectTable.AOSO_ID))

    % Assemble the corneal curvature, spherical error, and axial length.
    % For some subjects there are missing measures so we model the default
    % value for this eye.
    cc=[subjectTable.K1_average(ii),subjectTable.K2_average(ii),subjectTable.K1_angle_average(ii)];
    if any(isnan(cc))
        cc = [];
    end
    SR = subjectTable.Spherical_Error_average(ii);
    if isnan(SR)
        SR=[];
    end
    axialLength = subjectTable.Axial_Length_average(ii);
    if isnan(SR)
        axialLength=[];
    end

    % Create the model eye
    eye = modelEyeParameters('axialLength',axialLength,'sphericalAmetropia',SR,'measuredCornealCurvature',cc);
    
    % Extract the quadric surface for the vitreo-retinal interface
    S = eye.retina.S;
    
    % Define the visual field domain over which we will make the measure
    horizVals = -15:2.5:15;
    vertVals = -15:2.5:15;

    % Define an empty matrix to hold the results
    mmPerDeg = nan(length(horizVals),length(vertVals));

    % Define the delta deg
    deltaDegEuclidean = 1;
    deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2) 0];
    
    % Loop over horizontal and vertical field positions
    for jj = 1:length(horizVals)
        for kk = 1:length(vertVals)
            % The position in the field relative to the optical axis of the
            % eye
            degField = [horizVals(jj) vertVals(kk) 0] + p.Results.alpha;
            % Obtain the retinal points that are delta degrees on either
            % side of the specified degree field position
            [~,X0,angleError0] = findRetinaFieldPoint( eye, degField - deltaAngles./2);
            [~,X1,angleError1] = findRetinaFieldPoint( eye, degField + deltaAngles./2);
            % If the ray trace was accurate, calculate and store the
            % distance
            if angleError0 < 1e-3 && angleError1 < 1e-3
                % This is the minimal geodesic distance across the
                % ellipsoidal surface (like the great circle on a sphere).
                % Not using this because of too many failures across the
                % umbilical point, which is within the sampled area
                %{
                distance = abs(quadric.panouGeodesicDistance(S,[],[],X0,X1));
                %}
                % This is the Euclidean distance between the points
                distance = sqrt(sum((X0-X1).^2));
                % Need to divide by the delta distance to express as mm per
                % degree of visual angle
                mmPerDeg(jj,kk) = distance / deltaDegEuclidean;
            end
        end
    end
    % Give some console update
    fprintf(['Done subject ' num2str(subjectTable.AOSO_ID(ii)) '\n']);
    % Store the map in a cell array that accumulates across the parfor
    resultSet(ii) = {mmPerDeg};
end

% Write out the maps
for ii = 1:length(subjectTable.AOSO_ID)
    outfile = fullfile(saveDir,[num2str(subjectTable.AOSO_ID(ii)) '_mmPerDegMap.mat']);
    mmPerDeg = resultSet{ii};
    save(outfile,'mmPerDeg');
end
