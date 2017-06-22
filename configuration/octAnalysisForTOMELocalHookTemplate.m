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

%% Say hello
fprintf('Running OLFlickerSensitivitylocal hook\n');

%% Set preferences

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
    case 'connectome'
        dropboxBaseDir = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)'];
        dataPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/TOME_data/'];
    otherwise
        dropboxBaseDir = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)'];
        dataPath = ['/Users/' userID '/Dropbox (Aguirre-Brainard Lab)/MELA_data/'];
end

% Set the Dropox path
setpref('octAnalysisForTOME', 'dropboxPath', dropboxBaseDir);

% Set the data path
setpref('OneLight', 'dataPath', dataPath);

% Set the modulation path
setpref('OneLight', 'modulationPath', fullfile(dropboxBaseDir, 'MELA_materials', 'modulations/'));

% Set the materials path
setpref('OneLight', 'materialsPath', fullfile(dropboxBaseDir, 'MELA_materials/'));

% Set the cache path
setpref('OneLight', 'cachePath', fullfile(dropboxBaseDir, 'MELA_materials', 'cache/'));

% Set the cache path
setpref('OneLight', 'OneLightCalData', fullfile(dropboxBaseDir, 'MELA_materials', 'OneLightCalData/'));

% Set the default speak rate
setpref('OneLight', 'SpeakRateDefault', 230);

% Add OmniDriver.jar to java path
OneLightDriverPath = tbLocateToolbox('OneLightDriver');
JavaAddToPath(fullfile(OneLightDriverPath,'xOceanOpticsJava/OmniDriver.jar'),'OmniDriver.jar');

% Point at the code
olFlickerProjectDir = tbLocateProject('OLFlickerSensitivity');
setpref('OneLight', 'OLFlickerSensitivityBaseDir', fullfile(olFlickerProjectDir,'code'));

% Add OLFlickerSensitivity to the path
addpath(genpath(olFlickerProjectDir));

%% Experiment specific prefs
%
% MaxPulsePsychophysics
%
% set the nominal primaries path
setpref('MaxPulsePsychophysics','ModulationNominalPrimariesDir',fullfile(getpref('OneLight', 'materialsPath'),'Experiments','MaxPulsePsychophysics','ModulationNominalPrimaries'));
% set the spectrum sought primaries path
setpref('MaxPulsePsychophysics','ModulationCorrectedPrimariesDir',fullfile(getpref('OneLight', 'materialsPath'),'Experiments','MaxPulsePsychophysics','ModulationCorrectedPrimaries'));
% set the starts/stops path
setpref('MaxPulsePsychophysics','ModulationStartsStopsDir',fullfile(getpref('OneLight', 'materialsPath'),'Experiments','MaxPulsePsychophysics','ModulationStartsStops'));
% set the configuration files path
setpref('MaxPulsePsychophysics','ModulationConfigFilesDir',fullfile(getpref('OneLight', 'materialsPath'),'Experiments','MaxPulsePsychophysics','ModulationConfigFiles'));
