%% fitModelToAllData
% This routine loads the outputs of processDensityMaps.m and then fits the
% data with the polar cone density surface model.

% The overall result directory
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';

% Identify the aggregate data files
mergedFiles = dir([sourceDir '*_merged.mat']);
subNames = strrep(extractfield(mergedFiles,'name'),'_merged.mat','');

% Define some constants
supportLength = 1799;
imRdim = (supportLength+1)/2;
maxSupportDeg = 15;
supportDegDelta = 0.00773;

% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLength-1);


%% Loop through subjects and create the composite polar density image

dataMat = nan(supportLength,supportLength,length(subNames));
missingMerged = false(length(subNames));

for ss = 1:length(subNames)
    
    % A matrix to hold the data for this subject
    y = nan(supportLength,supportLength);
    
    % Load the aggregate data file
    merFile = fullfile(sourceDir,[subNames{ss} '_merged.mat']);
    if isfile(merFile)
        load(merFile,'data');
        y(:,:) = data.polarDensity(:,:);
    else
        str = fprintf(['No merged data for ' subNames{ss} ]);
        error(str);
    end        
    
    % Make sure that the supportDegDelta is as expected
    assert(abs(supportDegDelta-data.meta.supportDegDelta)<0.001);

    % Filter out any negative values
    y(y<0)=nan;
    dataMat(:,:,ss)=y;
    
end

% Fit the mean.
Y = nanmean(dataMat,3);
w = sum(~isnan(dataMat),3);
[p0, Yfit, fVal] = fitDensitySurface(Y,w,false,false,true,false);

%% Fit each subject with the reduced model
pSet = nan(20,length(subNames));
YfitSet = nan(size(dataMat));
YResidualSet = nan(size(dataMat));
fValSet= nan(1,length(subNames));
RSquaredSet= nan(4,length(subNames));
nonlconSet = nan(1,length(subNames));
polarThetaSet = nan(1,length(subNames));
polarMultiplierSet = nan(1,length(subNames));

fprintf('fitting...');
w1 = ones(size(Y));
for ii = 1:length(subNames)
    Y1 = squeeze(dataMat(:,:,ii));
    fprintf([num2str(ii),'...']);
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii), RSquaredSet(:,ii), nonlconSet(ii), polarThetaSet(ii), polarMultiplierSet(ii)] = fitDensitySurface(Y1,w1,true,true,false,false,p0);
    YResidualSet(:,:,ii) = Y1 - squeeze(YfitSet(:,:,ii));
end
fprintf('done\n');

% Save the individual subject fits
individualFitFile = fullfile(sourceDir,'individualSubjectFits.mat');
save(individualFitFile,'p0','pSet','YfitSet','fValSet','RSquaredSet','polarThetaSet','polarMultiplierSet','dataMat','subNames','YResidualSet')


%% Plot individual subject diagnostic plots
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];
polarRatio = (size(dataMat,1)+1)/360;

for ss=1:length(subNames)
    
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
        idx = find(supportDeg>=ii,1);
        if isempty(idx)
            density = coneDensityModel(ii,0,maxSupportDeg,pSet(:,ss));
            semilogy(repmat(density,1,size(dataMat,1)),'-r');
            text(300*polarRatio,density,sprintf('%2.1f°',ii));
        else
            semilogy(Y(:,idx),'.');
            hold on
            semilogy(Yfit(:,idx),'-r');
            text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
        end
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

