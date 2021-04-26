

% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')


% Load the confocal and split data
splitFiles=dir('densityAnalysis/*_split.mat');
splitNames = strrep(extractfield(splitFiles,'name'),'_split.mat','');

confocalFiles=dir('densityAnalysis/*_confocal.mat');
confocalNames = strrep(extractfield(confocalFiles,'name'),'_confocal.mat','');

supportLength = 1799;
maxSupportDeg = 15;
supportDegDelta = 0.0078;
conStart = 0.5;
conStop = 1.5;
splitStart = 1.75;

% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLength-1);
idxA = find(supportDeg>=conStart,1);
idxB = find(supportDeg>=conStop,1);
idxC = find(supportDeg>=splitStart,1);


%% Loop through subjects and create the composite polar density image

subNames = unique([splitNames confocalNames]);
dataMat = nan(supportLength,supportLength,length(subNames));

for ss = 1:length(subNames)
    
    % A matrix to hold the data for this subject
    y = nan(supportLength,supportLength);
    
    % Add the confocal data
    conFile = fullfile('densityAnalysis',[subNames{ss} '_confocal.mat']);
    if isfile(conFile)
        load(conFile,'data');
        y(:,idxA:idxB) = data.polarDensity(:,idxA:idxB);
    end
    
    % Add the split data
    splitFile = fullfile('densityAnalysis',[subNames{ss} '_split.mat']);
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
    end
    
    % Filter out any negative values
    y(y<0)=nan;
    dataMat(:,:,ss)=y;
    
end


% Fit the mean
Y = nanmean(dataMat,3);
w = sum(~isnan(dataMat),3);
[p0, Yfit, fVal] = fitDensitySurface(Y,w,false,false);


%% plot

% Values needed for the plots
polarRatio = (size(dataMat,1)+1)/360;
nFourier = 4;
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];

figure
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);
surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none');
hold on
plot3(X(:),P(:),Y(:),'.k')
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');


figure
for ii = [0.375 0.75 1.5 3 6 10]
    idx = find(supportDeg>ii,1);
    semilogy(Y(:,idx),'.');
    hold on
    semilogy(Yfit(:,idx),'-r');
    text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
end
xticks(meridianAngles*polarRatio);
xticklabels(meridianLabels);
ylim([10^2,10^4]);
ylabel('log_1_0 density [cones/deg^2]')


figure
for mm=1:4
    subplot(2,2,mm)
    plot(supportDeg,Y(round((meridianAngles(mm))*polarRatio+1),:),'.k');
    hold on
    plot(supportDeg,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
    xlabel('Eccentricity [deg]');
    ylabel('Density [cones/deg^2]');
    title(meridianLabels{mm});
end

figure
meridianSpec = {'-r','--b','-m','--g'};
for mm=1:4
    plot(supportDeg,Yfit(meridianAngles(mm)*polarRatio+1,:),meridianSpec{mm},'LineWidth',2);
    hold on
end
xlabel('Eccentricity [deg]');
ylabel('Density [cones/deg^2]');
legend(meridianLabels(1:4));


figure
for gg = 1:nFourier
    ph = p((gg-1)*4+5);
    f1 = p((gg-1)*4+6);
    f2 = p((gg-1)*4+7);
    f3 = p((gg-1)*4+8);
    g = f1.*gampdf(supportDeg,f2,f3)./max(gampdf(0:0.01:maxSupportDeg,f2,f3));
    plot(supportDeg,g)
    hold on
end
lengendLabels = {'cos1','sin1','cos2','cos4'};
legend(lengendLabels(1:nFourier));
xlabel('Eccentricity [deg]');
ylabel('Modulation [a.u.]');


figure
imagesc(w)
axis square
yticks(meridianAngles*polarRatio);
yticklabels(meridianLabels);
deltaDeg = supportDeg(2)-supportDeg(1);
xticks(round(0:1/deltaDeg:floor(maxSupportDeg)/deltaDeg));
xticklabels(0:1:floor(maxSupportDeg));
xlabel('Eccentricity [deg]');
colorbar
title('Weight map');



%% Fit each subject with the reduced model
pSet = nan(20,length(subNames));
YfitSet = nan(size(dataMat));
fValSet= nan(1,length(subNames));

fprintf('fitting...');
w1 = ones(size(Y));
for ii = 1:length(subNames)
    Y1 = squeeze(dataMat(:,:,ii));
    fprintf([num2str(ii),'...']);
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii)] = fitDensitySurface(Y1,w1,true,true,p0);
end
fprintf('done\n');


%% Plot individual subject fits
for ss=1:length(subNames)
    Y = squeeze(dataMat(:,:,ss));
    Yfit = squeeze(YfitSet(:,:,ss));
    p = pSet(:,ss);
        
    figure
    for ii = [0.375 0.75 1.5 3 6 10]
        idx = find(supportDeg>ii,1);
        semilogy(Y(:,idx),'.');
        hold on
        semilogy(Yfit(:,idx),'-r');
        text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
    end
    xticks(meridianAngles*polarRatio);
    xticklabels(meridianLabels);
    ylim([10^2,10^4]);
    ylabel('log_1_0 density [cones/deg^2]')
    
    
    figure
    for mm=1:4
        subplot(2,2,mm)
        plot(supportDeg,Y(round((meridianAngles(mm))*polarRatio+1),:),'.k');
        hold on
        plot(supportDeg,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
        xlabel('Eccentricity [deg]');
        ylabel('Density [cones/deg^2]');
        title(meridianLabels{mm});
    end
    
    pause
end

