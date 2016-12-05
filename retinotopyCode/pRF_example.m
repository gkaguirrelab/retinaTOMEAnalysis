%   This script provides example code to load data from an example TOME
%   subject, calculate population receptive fields, and visualize the
%   eccentricty and polar angle maps on the cortical surface.
%
%   Written by Andrew S Bock Oct 2016

%% Run on the UPenn cluster
% params.sessionDir   = '/data/jag/TOME/TOME_3005/100316';
% params.logDir       = '/data/jag/TOME/LOGS';
% makePRFshellScripts(params)

%% Run locally
% Set inputs
subjectName             = 'TOME_3005';
params.sessionDir       = '/Users/michael2/Documents/TOME/TOME_3005/100316/';
params.outDir           = fullfile(params.sessionDir,'pRFs');
hemis                   = {'lh','rh'};
% Get the retinotopy runs and stimulus files
b                       = listdir(fullfile(params.sessionDir,'*RETINO*'),'dirs');
stimFiles               = listdir(fullfile(params.sessionDir,'Stimuli','*RETINO*'),'files');
% run `makePRFmaps
for i = 2:length(b)
    for hh = 1:length(hemis)
        thisStim        = find(~cellfun('isempty',strfind(stimFiles,sprintf('run%02d',i))));
        params.stimFile = fullfile(params.sessionDir,'Stimuli',stimFiles{thisStim});
        params.inVol    = fullfile(params.sessionDir,b{i},['wdrf.tf.surf.' hemis{hh} '.nii.gz']);
        params.baseName         = sprintf([hemis{hh} '.run%02d'],i);
        % Calculate pRFs, save maps
        pRFs                    = makePRFmaps(params);
        % Plot the pRFs
        % Threshold by fit
        goodInd                 = pRFs.co>=sqrt(0.05);
        ecc                     = pRFs.ecc;
        pol                     = pRFs.pol;
        co                      = pRFs.co;
        sig                     = pRFs.sig;
        ecc(~goodInd)           = nan;
        pol(~goodInd)           = nan;
        co(~goodInd)            = nan;
        sig(~goodInd)           = nan;
        % Visualize maps
%         surface_plot('ecc',ecc,subjectName);
%         surface_plot('pol',pol,subjectName);
%         surface_plot('co',co,subjectName);
%         surface_plot('sig',sig,subjectName);
    end
end
