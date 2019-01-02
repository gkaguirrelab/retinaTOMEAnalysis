% Convert OCT seg XML to .mat volumes


subjectsToProcess={'11083','11084','11086','11087','11088','11089','11091','11092','11093','11095','11096','11097','11098','11099'};

% set dropbox directory
[~,hostname] = system('hostname');
hostname = strtrim(lower(hostname));

% handle hosts with custom dropbox locations
switch hostname
    case 'seele.psych.upenn.edu'
        dropboxDir = '/Volumes/seeleExternalDrive/Dropbox (Aguirre-Brainard Lab)';
    case 'magi-1-melchior.psych.upenn.edu'
        dropboxDir = '/Volumes/melchiorExternalDrive/Dropbox (Aguirre-Brainard Lab)';
    case 'magi-2-balthasar.psych.upenn.edu'
        dropboxDir = '/Volumes/balthasarExternalDrive/Dropbox (Aguirre-Brainard Lab)';
    otherwise
        [~, userName] = system('whoami');
        userName = strtrim(userName);
        dropboxDir = ...
            fullfile('/Users', userName, ...
            'Dropbox (Aguirre-Brainard Lab)');
end


dataRootDir = [dropboxDir '/AOSO_analysis/OCTExplorerSegmentationData'];
inFileSuffix = '_Surfaces_Retina-JEI-Final.xml';
outFileSuffix = '_Surfaces_Retina-JEI-Final.mat';

% Obtain the paths to all of the pupil data files within the specified
% directory, including within sub-drectories.
fileListStruct=subdir(fullfile(dataRootDir,['*' inFileSuffix]));

% If we found at least one pupil data file, then proceed.
if ~isempty(fileListStruct)
    
    for ii = 1:length(fileListStruct)
        
        processThisOneFlag = true;
        % Check if we have a list of subjects to process
        if ~isempty(subjectsToProcess)
            % Check if this subject is on the list
            if ~contains(fileListStruct(ii).name,subjectsToProcess)
                processThisOneFlag = false;
            end
        end
        
        if processThisOneFlag
            % Convert this XML file to a volume
            mask = xml2volmask(fileListStruct(ii).name);
            
            % Create an output file path
            outFileRoot = extractBefore(fileListStruct(ii).name,inFileSuffix);
            outputFile = [outFileRoot outFileSuffix];
            
            % Save the file
            save(outputFile,'mask');
            
            % Report the conversion
            fprintf([outputFile '\n']);
        end
    end
end




function varargout = subdir(varargin)
%SUBDIR Performs a recursive file search
%
% subdir
% subdir(name)
% files = subdir(...)
%
% This function performs a recursive file search.  The input and output
% format is identical to the dir function.
%
% Input variables:
%
%   name:   pathname or filename for search, can be absolute or relative
%           and wildcards (*) are allowed.  If ommitted, the files in the
%           current working directory and its child folders are returned
%
% Output variables:
%
%   files:  m x 1 structure with the following fields:
%           name:   full filename
%           date:   modification date timestamp
%           bytes:  number of bytes allocated to the file
%           isdir:  1 if name is a directory; 0 if no
%
% Example:
%
%   >> a = subdir(fullfile(matlabroot, 'toolbox', 'matlab', '*.mat'))
%
%   a =
%
%   67x1 struct array with fields:
%       name
%       date
%       bytes
%       isdir
%
%   >> a(2)
%
%   ans =
%
%        name: '/Applications/MATLAB73/toolbox/matlab/audiovideo/chirp.mat'
%        date: '14-Mar-2004 07:31:48'
%       bytes: 25276
%       isdir: 0
%
% See also:
%
%   dir

% Copyright 2006 Kelly Kearney


%---------------------------
% Get folder and filter
%---------------------------

narginchk(0,1);
nargoutchk(0,1);

if nargin == 0
    folder = pwd;
    filter = '*';
else
    [folder, name, ext] = fileparts(varargin{1});
    if isempty(folder)
        folder = pwd;
    end
    if isempty(ext)
        if isdir(fullfile(folder, name))
            folder = fullfile(folder, name);
            filter = '*';
        else
            filter = [name ext];
        end
    else
        filter = [name ext];
    end
    if ~isdir(folder)
        error('Folder (%s) not found', folder);
    end
end

%---------------------------
% Search all folders
%---------------------------

pathstr = genpath_local(folder);
pathfolders = regexp(pathstr, pathsep, 'split');  % Same as strsplit without the error checking
pathfolders = pathfolders(~cellfun('isempty', pathfolders));  % Remove any empty cells

Files = [];
pathandfilt = fullfile(pathfolders, filter);
for ifolder = 1:length(pathandfilt)
    NewFiles = dir(pathandfilt{ifolder});
    if ~isempty(NewFiles)
        fullnames = cellfun(@(a) fullfile(pathfolders{ifolder}, a), {NewFiles.name}, 'UniformOutput', false);
        [NewFiles.name] = deal(fullnames{:});
        Files = [Files; NewFiles];
    end
end

%---------------------------
% Prune . and ..
%---------------------------

if ~isempty(Files)
    [~, ~, tail] = cellfun(@fileparts, {Files(:).name}, 'UniformOutput', false);
    dottest = cellfun(@(x) isempty(regexp(x, '\.+(\w+$)', 'once')), tail);
    Files(dottest & [Files(:).isdir]) = [];
end

%---------------------------
% Output
%---------------------------

if nargout == 0
    if ~isempty(Files)
        fprintf('\n');
        fprintf('%s\n', Files.name);
        fprintf('\n');
    end
elseif nargout == 1
    varargout{1} = Files;
end

end % subdir

function [p] = genpath_local(d)
% Modified genpath that doesn't ignore:
%     - Folders named 'private'
%     - MATLAB class folders (folder name starts with '@')
%     - MATLAB package folders (folder name starts with '+')

files = dir(d);
if isempty(files)
    return
end
p = '';  % Initialize output

% Add d to the path even if it is empty.
p = [p d pathsep];

% Set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));
dirs = files(isdir);  % Select only directory entries from the current listing

for i=1:length(dirs)
    dirname = dirs(i).name;
    if    ~strcmp( dirname,'.') && ~strcmp( dirname,'..')
        p = [p genpath(fullfile(d,dirname))];  % Recursive calling of this function.
    end
end

end % genpath_local
