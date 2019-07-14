function retinaTOMEAnalysisLocalHook
% retinaTOMEAnalysisLocalHook - Configure things for working on OneLight projects.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'retinaTOMEAnalysisConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
projectName = 'retinaTOMEAnalysis';
 
%% Say hello
fprintf('Running % local hook\n',projectName);
 
%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end
 
%% Specify project location
projectBaseDir = tbLocateProject(projectName);

% Obtain the Dropbox path
[~,hostname] = system('hostname');
hostname = strtrim(lower(hostname));

% handle hosts with custom dropbox locations
switch hostname
    case 'seele.psych.upenn.edu'
        dropboxBaseDir = '/Volumes/seeleExternalDrive/Dropbox (Aguirre-Brainard Lab)';
    case 'magi-1-melchior.psych.upenn.edu'
        dropboxBaseDir = '/Volumes/melchiorBayTwo/Dropbox (Aguirre-Brainard Lab)';
    case 'magi-2-balthasar.psych.upenn.edu'
        dropboxBaseDir = '/Volumes/balthasarExternalDrive/Dropbox (Aguirre-Brainard Lab)';
    otherwise
        [~, userName] = system('whoami');
        userName = strtrim(userName);
        dropboxBaseDir = ...
            fullfile('/Users', userName, ...
            'Dropbox (Aguirre-Brainard Lab)');
end

%% Set preferences for project output
setpref(projectName,'dropboxBaseDir',dropboxBaseDir); % main directory path 
setpref(projectName,'projectBaseDir',projectBaseDir); % main directory path 

%% Set the vlfeatroot
% vlfeat is a compiled library of vision algorithms. Download the latest
% version from vlfeat.org. Unpack within userpath/src/. Update the version
% ID in the vlfeatroot variable below.
VLFEATROOT = fullfile(userpath(),'src','vlfeat-0.9.21');
run(fullfile(VLFEATROOT,'toolbox','vl_setup'));
fprintf('vlfeat version: ');
vl_version

end
