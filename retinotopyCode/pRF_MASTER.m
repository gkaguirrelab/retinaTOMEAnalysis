% MASTER script to process the pRF data for the octAnalysisForTOME project


%% Set various inputs
sessions = {...
    '/Users/michael2/Documents/TOME/TOME_3003/091616/', ...
    '/Users/michael2/Documents/TOME/TOME_3005/100316/' ...
    };
dbDir = '/Users/michael2/Dropbox-Aguirre-Brainard-Lab/';
dbOutDirs = {...
    fullfile(dbDir,'TOME_3003','091616','pRFs'),...
    fullfile(dbDir,'TOME_3005','100316','pRFs') ...
    };
hemis                   = {'lh','rh'};
%% Create pRF maps
for i = 1:length(sessions)
    params.sessionDir       = sessions{i};
    params.outDir           = fullfile(params.sessionDir,'pRFs');
    % Get the retinotopy runs and stimulus files
    b                       = listdir(fullfile(params.sessionDir,'*RETINO*'),'dirs');
    stimFiles               = listdir(fullfile(params.sessionDir,'Stimuli','*RETINO*'),'files');
    % run `makePRFmaps
    for j = 1:length(b)
        thisStim            = find(~cellfun('isempty',strfind(stimFiles,sprintf('run%02d',j))));
        for k = 1:length(hemis)
            params.stimFile = fullfile(params.sessionDir,'Stimuli',stimFiles{thisStim});
            params.inVol    = fullfile(params.sessionDir,b{j},['wdrf.tf.surf.' hemis{k} '.nii.gz']);
            params.baseName = sprintf([hemis{k} '.run%02d'],j);
            % Calculate pRFs, save maps
            pRFs            = makePRFmaps(params);
        end
    end
end
%% Average the maps
for i = 1:length(sessions)
    params.inDir            = fullfile(sessions{i},'pRFs');
    params.outDir           = params.inDir;
    for j = 1:length(hemis)
        params.baseName     = hemis{j};
        avgPRFmaps(params)
    end
end
%% Copy files to Dropbox
for i = 1:length(sessions)
    inDir            = fullfile(sessions{i},'pRFs');
    outDir           = dbOutDirs{i};
    system(['cp ' fullfile(inDir,'*.nii.gz') ' ' outDir]);
end