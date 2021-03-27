function makeModelEyes(saveDir, varargin)
% Creates and saves a model eye for each subject
%
% Description:
%   This will subsequently be used to convert between degrees and mm on the
%   retina.
%
% Examples:
%{
    dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
    saveDir = fullfile(dropboxBaseDir,'AOSO_analysis','eyeModels');
    makeModelEyes(saveDir)
%}


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

eyeModels = cell(1,length(subjectTable.AOSO_ID));

for ii = 1:length(subjectTable.AOSO_ID)

    % Assemble the corneal curvature, spherical error, and axial length.
    axialLength = subjectTable.Axial_Length_average(ii);
    SR = subjectTable.Spherical_Error_average(ii);
    kvals = subjectTable.kvals{ii};
    if ~isempty(kvals)
        kvals = eval(kvals);
    else
        kvals = [];
    end

    % Create the model eye
    eyeModels{ii} = modelEyeParameters('axialLength',axialLength,'sphericalAmetropia',SR,'kvals',kvals);
    
    % Give some console update
    fprintf(['Done subject ' num2str(subjectTable.AOSO_ID(ii)) '\n']);
end

% Write out the results
outfile = fullfile(saveDir,'eyeModels.mat');
save(outfile,'eyeModels');

end % Main


