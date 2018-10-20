function OLFlickerSensitivityLocalHook
% OLFlickerSensitivityLocalHook - Configure things for working on OneLight projects.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'OLFlickerSensitivityConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
projectName = 'octAnalysisForTOME';
 
%% Say hello
fprintf('Running % local hook\n',projectName);
 
%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end
 
%% Specify project location
projectBaseDir = tbLocateProject(projectName);

% Obtain the Dropbox path
[~, userID] = system('whoami');
userID = strtrim(userID);
switch userID
    case {'melanopsin' 'pupillab'}
        dropboxBaseDir = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/'];
        dataPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
    case {'dhb'}
        dropboxBaseDir = ['/Users1'  '/Dropbox (Aguirre-Brainard Lab)/'];
        dataPath = ['/Users1/' '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];        
    otherwise
        dropboxBaseDir = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)'];
        DropBoxDataPath = [dropboxBaseDir '/retData/'];
        LocalDataPath = [projectBaseDir '/data'];
end
 
%% Set preferences for project output

setpref(projectName,'mainDir',projectBaseDir); % main directory path 
setpref(projectName,'DropBoxDataPath',DropBoxDataPath); % path to data stroed on dropbox
setpref(projectName,'LocalDataPath',LocalDataPath); % path to small file within the git repo 
end
