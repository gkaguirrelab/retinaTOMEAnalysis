%% fitModelToAllData
% This routine loads the outputs of processDensityMaps.m, generates
% composites of the data (confocal, split, and "fovea") and then fits the
% data with the polar cone density surface model.

% The overal result directory
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';

% Identify the confocal, split, and "fovea" data files
splitFiles=dir([sourceDir '*_split.mat']);
splitNames = strrep(extractfield(splitFiles,'name'),'_split.mat','');

confocalFiles=dir([sourceDir '*_confocal.mat']);
confocalNames = strrep(extractfield(confocalFiles,'name'),'_confocal.mat','');

foveaFiles=dir([sourceDir '*_fovea.mat']);
foveaNames = strrep(extractfield(confocalFiles,'name'),'_fovea.mat','');

% Define some constants
supportLength = 1799;
imRdim = (supportLength+1)/2;
maxSupportDeg = 15;
supportDegDelta = 0.0078;
conStart = 0.75;
conStop = 1.5;
splitStart = 1.75;
foveaStart = 0.3;
foveaStop = 0.5;

% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLength-1);
idxA = find(supportDeg>=conStart,1);
idxB = find(supportDeg>=conStop,1);
idxC = find(supportDeg>=splitStart,1);
idxD = find(supportDeg>=foveaStart,1);
idxE = find(supportDeg>=foveaStop,1);


%% Loop through subjects and create the composite polar density image

subNames = unique([splitNames confocalNames]);
dataMat = nan(supportLength,supportLength,length(subNames));
missingSplit = false(length(subNames));
missingFovea = false(length(subNames));

for ss = 1:length(subNames)
    
    % A matrix to hold the data for this subject
    y = nan(supportLength,supportLength);
    
    % Add the confocal data
    conFile = fullfile(sourceDir,[subNames{ss} '_confocal.mat']);
    if isfile(conFile)
        load(conFile,'data');
        y(:,idxA:idxB) = data.polarDensity(:,idxA:idxB);
    end
    
    % Add the split data
    splitFile = fullfile(sourceDir,[subNames{ss} '_split.mat']);
    hasSplit = false;
    if isfile(splitFile)
        load(splitFile,'data');
        splitBit = data.polarDensity(:,idxC:end);
        y(:,idxC:end) = splitBit;
        if sum(~isnan(splitBit(:)))>0
            hasSplit = true;
        end
    end
    if ~hasSplit
        fprintf(['No split data for ' subNames{ss} '\n']);
        missingSplit(ss) = true;
    end
    
    % Add the fovea data
    foveaFile = fullfile(sourceDir,[subNames{ss} '_fovea.mat']);
    hasFovea = false;
    if isfile(foveaFile)
        load(foveaFile,'data');
        foveaBit = data.polarDensity(:,1:idxE);
        y(:,1:idxE) = foveaBit;
        if sum(~isnan(foveaBit(:)))>0
            hasFovea = true;
        end
    end    
    if ~hasFovea
        fprintf(['No fovea data for ' subNames{ss} '\n']);
        missingFovea(ss) = true;
    end
    
    % Filter out any negative values
    y(y<0)=nan;
    dataMat(:,:,ss)=y;
    
end


% Fit the mean. We remove value very close to the fovea as not all subjects
% have these, leading to a bias of the mean close to the fovea towards
% those subjects with lower cone densities.
Y = nanmean(dataMat,3);
w = sum(~isnan(dataMat),3);
Y(:,1:idxD)=nan;
w(:,1:idxD)=0;
[p0, Yfit, fVal] = fitDensitySurface(Y,w,false,false,true);



%% Fit each subject with the reduced model
pSet = nan(20,length(subNames));
YfitSet = nan(size(dataMat));
fValSet= nan(1,length(subNames));

fprintf('fitting...');
w1 = ones(size(Y));
for ii = 1:length(subNames)
    if missingSplit(ii)
        continue
    end
    Y1 = squeeze(dataMat(:,:,ii));
    fprintf([num2str(ii),'...']);
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii), polarMultiplierSet(ii)] = fitDensitySurface(Y1,w1,true,true,false,p0);
end
fprintf('done\n');

% Save the individual subject fits
individualFitFile = fullfile(sourceDir,'individualSubjectFits.mat');
save(individualFitFile,'pSet','YfitSet','fValSet','polarMultiplierSet','dataMat','subNames')


%% Plot individual subject diagnostic plots
for ss=1:length(subNames)
    
    if missingSplit(ii)
        continue
    end
    
    Y = squeeze(dataMat(:,:,ss));
    Yfit = squeeze(YfitSet(:,:,ss));
    p = pSet(:,ss);
    
    % Mean polar data and model fit
    figHandle = figure();
    X = repmat(supportDeg,supportLength,1);
    P = repmat(linspace(0,360,supportLength)',1,supportLength);
    surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none');
    hold on
    plot3(X(:),P(:),Y(:),'.k')
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [deg]');
    zlabel('Density [cones/deg^2]');
    view(45,15)
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_polarSurface.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);
        
    % Polar bands and model fit
    figHandle = figure();    
    for ii = [0.1875 0.375 0.75 1.5 3 6 10]
        idx = find(supportDeg>ii,1);
        semilogy(Y(:,idx),'.');
        hold on
        semilogy(Yfit(:,idx),'-r');
        text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
    end
    xticks(meridianAngles*polarRatio);
    xticklabels(meridianLabels);
    ylim([10^2,2*10^4]);
    ylabel('log_1_0 density [cones/deg^2]')
    
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_polarModelFit.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);
    
    % Plot by meridian
    figHandle = figure();
    
    for mm=1:4
        subplot(2,2,mm)
        plot(supportDeg,Y(round((meridianAngles(mm))*polarRatio+1),:),'.k');
        hold on
        plot(supportDeg,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
        xlabel('Eccentricity [deg]');
        ylabel('Density [cones/deg^2]');
        title(meridianLabels{mm});
    end
    
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_MeridianModelFit.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);    
    
end

