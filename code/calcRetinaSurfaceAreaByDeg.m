
% Load the subject data table
subjectTableFileName='/Volumes/balthasarExternalDrive/Dropbox (Aguirre-Brainard Lab)/TOME_subject/TOME-AOSO_SubjectInfo.xlsx';
opts = detectImportOptions(subjectTableFileName);
subjectTable = readtable(subjectTableFileName, opts);

% Set the save directory and initialize the cell array resultSet
saveDir = '/Volumes/balthasarExternalDrive/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis';
resultSet = {};

parfor (ii = 1:length(subjectTable.AOSO_ID))
%for ii = 1:length(subjectTable.AOSO_ID)

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
    
    % Set the alpha angles to the population mean
    alpha = [5.45 2.5 0];
    
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
            degField = [horizVals(jj) vertVals(kk) 0] + alpha;
            % Obtain the retinal points that are delta degrees on either
            % side of the specified degree field position
            [~,X0,angleError0] = findRetinaFieldPoint( eye, degField - deltaAngles./2);
            [~,X1,angleError1] = findRetinaFieldPoint( eye, degField + deltaAngles./2);
            % If the ray trace was accurate, calculate and store the
            % distance
            if angleError0 < 1e-3 && angleError1 < 1e-3
                % This is the minimal geodesic distance across the
                % ellipsoidal surface (like the great circle on a sphere)
                distance = abs(quadric.panouGeodesicDistance(S,[],[],X0,X1));
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
    outfile = fullfile(saveDir,'mmPerDegMaps',[num2str(subjectTable.AOSO_ID(ii)) '_mmPerDegMap.mat']);
    mmPerDeg = resultSet{ii};
    save(outfile,'mmPerDeg');
end
