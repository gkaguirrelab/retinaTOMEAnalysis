function makeLeftRightMontage(dataRootDir, saveDir, varargin)



%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('dataRootDir',@ischar);
p.addRequired('saveDir',@ischar);

% Optional analysis params
p.addParameter('layerSetLabels',{'RGCIPL','RNFL','OPL','TotalRetina'},@iscell);
p.addParameter('layerSets',[{2:3},{1},{6},{1:11}],@iscell);
p.addParameter('foveaThickThresh',75,@isscalar);
p.addParameter('showPlots',true,@islogical);
p.addParameter('verbose',true,@islogical);
p.addParameter('artifactTrim',[1280, 1536, 1280, 1536],@isscalar);


%% Parse and check the parameters
p.parse(dataRootDir, saveDir, varargin{:});

% Obtain a list of subjects
rawSubjectList = dir(dataRootDir);

% How many sets are there
nSets = length(p.Results.layerSets);

% Prepare a variable to hold the overlap maps
overlapMaps = struct;

% Loop over subjects
for ss=1:length(rawSubjectList)
    
    % Get the name of this subject
    subjectName = rawSubjectList(ss).name;
    
    % Assemble the path to the data directory for this subject
    dirPath = fullfile(rawSubjectList(ss).folder,subjectName);
    
    % Obtain the list of result files for this subject
    fileList = dir(fullfile(dirPath,'*.mat'));
    
    % If there are no files for this subject, move on to the next subject
    if isempty(fileList)
        continue
    end
    
    % If there are more than two result files for this subject, error.
    if length(fileList)>2
        error('Why are there more than two eyes for this subject?')
    end
    
    % Report this subject's name
    if p.Results.verbose
        fprintf(['Processing subject: ' subjectName '\n']);
    end
    
    % Open a figure to display the maps
    if p.Results.showPlots
        figure('Name',subjectName,'NumberTitle','off');
    end
    
    if length(fileList)==1
        warning(['subject ' subjectName ' has just one eye']);
    end
    
    % Loop over the right and left eye
    for ii = 1:length(fileList)
        
        % Assemble the name of this file and load it
        fileName = fullfile(fileList(ii).folder,fileList(ii).name);
        load(fileName,'subject');
        
        % Empty the variable that will hold the fovea coordinates
        foveaCoord=[];
        
        % Loop over the layer sets
        for jj=1:nSets
            
            % Get the thickness by summing the layers defined in
            % layerIdx
            thisThick = sum(subject.BothMeanLayerThicknessesOnSLOInterp(:,:,p.Results.layerSets{jj}),3);
            
            % Get the imageSize
            imageSize = size(thisThick);
            
            % If this is the first file and the first layer set, define a
            % variable to hold results across layer sets and eyes
            if ii==1 && jj==1
                everyThicknessMap = nan(nSets,2,imageSize(1),imageSize(2));
            end
            
            % If not yet defined, set up the overlapMaps variable for this
            % layer
            if ~isfield(overlapMaps,p.Results.layerSetLabels{jj})
                overlapMaps.(p.Results.layerSetLabels{jj}) = ...
                    ones(imageSize(1),imageSize(2));
            end
            
            % If these are data from the left eye, mirror reverse
            if contains(fileList(ii).name,'_OS.mat')
                thisThick = fliplr(thisThick);
            end
            
            % Set points with zero thickness to nan
            thisThick(thisThick==0)=nan;
            
            % If this is the first layer set, then we are working with the
            % RGC-IPL layer. Use this to find the fovea at the center of the
            % image. This is defined by taking the weighted mean of the
            % thinnest portion of the central 20% of the image
            if jj==1
                w=thisThick;
                w(1:round(imageSize(1)*0.4),:)=nan;
                w(round(imageSize(1)*0.6):end,:)=nan;
                w(:,1:round(imageSize(1)*0.4))=nan;
                w(:,round(imageSize(1)*0.6):end)=nan;
                w(w>p.Results.foveaThickThresh)=nan;
                w=max(max(w))-w;
                w=w(1:prod(imageSize));
                w=w./nansum(w);
                [X,Y] = ind2sub(imageSize,1:prod(imageSize));
                foveaCoord(1) = nansum(X.*w(1:prod(imageSize)));
                foveaCoord(2) = nansum(Y.*w(1:prod(imageSize)));
            end
            
            % Save a mask of the nan values
            nanMask = zeros(size(thisThick));
            nanMask(isnan(thisThick))=1;
            thisThick(isnan(thisThick))=0;
            
            % Shift the thickness and nanmask images to align the foveas
            thisThick = imtranslate(thisThick,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
            nanMask = imtranslate(nanMask,fliplr((imageSize./2)-foveaCoord),'method','cubic','FillValues',nan);
            
            % Restore the nans
            thisThick(nanMask>0.25)=nan;
            thisThick(thisThick==0)=nan;
            
            % Trim away artifacts in the area of the optic disc
            thisThick(:,p.Results.artifactTrim(jj):end)=nan;
            
            % Store the this thickness map in the full array
            if contains(fileList(ii).name,'_OS.mat')
                everyThicknessMap(jj,1,:,:)=thisThick;
            else
                everyThicknessMap(jj,2,:,:)=thisThick;
            end
            
            % Show the RGC+IPL map and the RNFL map
            if p.Results.showPlots
                switch p.Results.layerSetLabels{jj}
                    case 'OPL'
                        tmp = fliplr(thisThick);
                        if contains(fileList(ii).name,'_OS.mat')
                            subplot(2,2,1)
                            tmp = fliplr(tmp);
                        else
                            subplot(2,2,2)
                        end
                        imagesc(tmp);
                        hold on
                        plot(imageSize(1)/2,imageSize(2)/2,'+k')
                        axis square
                        drawnow
                    case 'TotalRetina'
                        tmp = fliplr(thisThick);
                        if contains(fileList(ii).name,'_OS.mat')
                            subplot(2,2,3)
                            tmp = fliplr(tmp);
                        else
                            subplot(2,2,4)
                        end
                        imagesc(tmp);
                        hold on
                        plot(imageSize(1)/2,imageSize(2)/2,'+k')
                        axis square
                        drawnow
                end
            end
            
        end % Loop over layer sets
        
    end % Looping over the two eyes
    
    % Obtain the nanmean across eyes
    everyThicknessMap = squeeze(nanmean(everyThicknessMap,2));
    
    % Loop over sets and assemble the data into a struct
    averageMaps = struct();
    for jj=1:nSets
        
        % Store the average map
        averageMaps.(p.Results.layerSetLabels{jj}) = ...
            squeeze(everyThicknessMap(jj,:,:));
        
        % Update the overlapMaps to nan those points that are nan here
        overlapMaps.(p.Results.layerSetLabels{jj})(isnan(averageMaps.(p.Results.layerSetLabels{jj}))) = nan;
    end
    
    % Save the structure variable
    outfileDir = fullfile(saveDir,subjectName);
    if ~exist(outfileDir,'dir')
        mkdir(outfileDir)
    end
    outfile = fullfile(outfileDir,[subjectName,'_averageMaps.mat']);
    save(outfile,'averageMaps');
    
end % Looping over subjects

% Save the overlapMaps
outfile = fullfile(saveDir,'overlapMaps.mat');
save(outfile,'overlapMaps');



