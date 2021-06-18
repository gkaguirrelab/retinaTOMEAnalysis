

% The overal result directory

sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';

% Load the confocal and split data
splitFiles=dir([sourceDir '*_split.mat']);
splitNames = strrep(extractfield(splitFiles,'name'),'_split.mat','');

confocalFiles=dir([sourceDir '*_confocal.mat']);
confocalNames = strrep(extractfield(confocalFiles,'name'),'_confocal.mat','');

foveaFiles=dir([sourceDir '*_fovea.mat']);
foveaNames = strrep(extractfield(confocalFiles,'name'),'_fovea.mat','');


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
[p0, Yfit, fVal] = fitDensitySurface(Y,w,false,false);


%% plot

% Values needed for the plots
polarRatio = (size(dataMat,1)+1)/360;
nFourier = 4;
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];


% Mean polar data and model fit
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
view(45,15)
plotFileName = fullfile(sourceDir,'figures','Fig01_meanModelFitPolar.pdf');
saveas(gcf,plotFileName);

% Illustrate the model surface components
figure
for cc = 1:4
    subplot(2,2,cc)
    switch cc
        case 1
            Ymodel = cosd(P+p0(5));
        case 2
            Ymodel = sind(P+p0(5));
        case 3
            Ymodel = cosd(P.*2+p0(13));
        case 4
            Ymodel = cosd(P.*4+p0(17));
    end
    surf(X,P,Ymodel,'FaceAlpha',0.5,'EdgeColor','none');
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [deg]');
    zlabel('Density modulation');
    view(45,15)
end
plotFileName = fullfile(sourceDir,'figures','Fig0X_polarModulationModel.pdf');
saveas(gcf,plotFileName);

% Illustrate the model surface components
figure
for cc = 1:4
    subplot(2,2,cc)
    switch cc
        case 1
            g = p0(6).*gampdf(X,p0(7),p0(8))./max(gampdf(0:0.01:maxSupportDeg,p0(7),p0(8)));
            Ymodel = g.*cosd(P+p0(5));
        case 2
            g = p0(10).*gampdf(X,p0(11),p0(12))./max(gampdf(0:0.01:maxSupportDeg,p0(11),p0(12)));
            Ymodel = g.*sind(P+p0(5));
        case 3
            g = p0(14).*gampdf(X,p0(15),p0(16))./max(gampdf(0:0.01:maxSupportDeg,p0(15),p0(16)));
            Ymodel = g.*cosd(P.*2+p0(13));
        case 4
            g = p0(18).*gampdf(X,p0(19),p0(20))./max(gampdf(0:0.01:maxSupportDeg,p0(19),p0(20)));
            Ymodel = g.*cosd(P.*4+p0(17));
    end
    surf(X,P,Ymodel,'FaceAlpha',0.5,'EdgeColor','none');
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [deg]');
    zlabel('Density modulation');
    view(45,15)
end
plotFileName = fullfile(sourceDir,'figures','Fig0X_polarModulationModelWithGamma.pdf');
saveas(gcf,plotFileName);

% Illustrate the model surface components
figure

g = p0(6).*gampdf(X,p0(7),p0(8))./max(gampdf(0:0.01:maxSupportDeg,p0(7),p0(8)));
Ymodel = g.*cosd(P+p0(5));
g = p0(10).*gampdf(X,p0(11),p0(12))./max(gampdf(0:0.01:maxSupportDeg,p0(11),p0(12)));
Ymodel = Ymodel+g.*sind(P+p0(5));
g = p0(14).*gampdf(X,p0(15),p0(16))./max(gampdf(0:0.01:maxSupportDeg,p0(15),p0(16)));
Ymodel = Ymodel+g.*cosd(P.*2+p0(13));
g = p0(18).*gampdf(X,p0(19),p0(20))./max(gampdf(0:0.01:maxSupportDeg,p0(19),p0(20)));
Ymodel = Ymodel+g.*cosd(P.*4+p0(17));
surf(X,P,Ymodel,'FaceAlpha',0.5,'EdgeColor','none');
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density modulation');
view(45,15)
plotFileName = fullfile(sourceDir,'figures','Fig0X_entirePolarModulation.pdf');
saveas(gcf,plotFileName);


% Mean Cartesian data and model fit
figure
cartSupportDeg=linspace(-max(supportDeg),max(supportDeg),imRdim);
[cartXDeg,cartYDeg] = meshgrid(cartSupportDeg,cartSupportDeg);
cartYfit = convertPolarMapToImageMap(Yfit,'imRdim',imRdim);
cartYfit(cartYfit<min([Y(:); Yfit(:)]))=nan;
cartY = convertPolarMapToImageMap(Y,'imRdim',imRdim);
cartY(cartY<min([Y(:); Yfit(:)]))=nan;
surf(cartXDeg,cartYDeg,cartYfit,'FaceAlpha',0.5,'EdgeColor','none');
hold on
%plot3(cartXDeg(:),cartYDeg(:),cartY(:),'.k')
xlabel('Eccentricity [deg]');
ylabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');
view(-140,22)
plotFileName = fullfile(sourceDir,'figures','Fig02_meanModelFitCartesian.pdf');
saveas(gcf,plotFileName);


% Polar angle density variation
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
plotFileName = fullfile(sourceDir,'figures','Fig03_densityVariationByPolarAngle.pdf');
saveas(gcf,plotFileName);


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
plotFileName = fullfile(sourceDir,'figures','Fig04_meanModelFitByMeridian.pdf');
saveas(gcf,plotFileName);


figure
meridianSpec = {'-r','--b','-m','--g'};
for mm=1:4
    plot(supportDeg,Yfit(meridianAngles(mm)*polarRatio+1,:),meridianSpec{mm},'LineWidth',2);
    hold on
end
xlabel('Eccentricity [deg]');
ylabel('Density [cones/deg^2]');
legend(meridianLabels(1:4));
plotFileName = fullfile(sourceDir,'figures','Fig05_modeledDensityByMeridian.pdf');
saveas(gcf,plotFileName);


figure
for gg = 1:nFourier
    ph = p0((gg-1)*4+5);
    f1 = p0((gg-1)*4+6);
    f2 = p0((gg-1)*4+7);
    f3 = p0((gg-1)*4+8);
    g = f1.*gampdf(supportDeg,f2,f3)./max(gampdf(0:0.01:maxSupportDeg,f2,f3));
    plot(supportDeg,g)
    hold on
end
lengendLabels = {'cos1','sin1','cos2','cos4'};
legend(lengendLabels(1:nFourier));
xlabel('Eccentricity [deg]');
ylabel('Modulation [proportion]');
plotFileName = fullfile(sourceDir,'figures','Fig06_gammaPDFWeightFunctions.pdf');
saveas(gcf,plotFileName);



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
plotFileName = fullfile(sourceDir,'figures','Fig07_nSubjectWeightMapPolar.pdf');
saveas(gcf,plotFileName);



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
    [pSet(:,ii), YfitSet(:,:,ii), fValSet(ii)] = fitDensitySurface(Y1,w1,true,true,p0);
end
fprintf('done\n');

% Save the individual subject fits
individualFitFile = fullfile(sourceDir,'individualSubjectFits.mat');
save(individualFitFile,'pSet','YfitSet','fValSet')



%% Plot the Curcio values
figure
figure
meridianSpec = {'-r','--b','-m','--g'};
for mm=1:4
    CurcioFitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(meridianAngles(mm));
    plot(supportDeg,CurcioFitConeDensitySqDegVisual(supportDeg),meridianSpec{mm});
    hold on
end
xlabel('Eccentricity [deg]');
ylabel('Density [cones/deg^2]');
legend(meridianLabels(1:4));
plotFileName = fullfile(sourceDir,'figures','Fig0x_CurcioValues.pdf');
saveas(gcf,plotFileName);



%% Plot the fit on the nasal meridian for each subject
figure
for meridianIdx = 1:4
    subplot(2,2,meridianIdx)
    % Add the Curcio profile along the specified meridian
    CurcioFitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(meridianAngles(meridianIdx));
    loglog(supportDeg,CurcioFitConeDensitySqDegVisual(supportDeg),'-r','LineWidth',5);
    hold on
    for ss = 1:length(subNames)
        if missingSplit(ss)
            continue
        end
        Yfit = squeeze(YfitSet(round((meridianAngles(meridianIdx))*polarRatio+1),:,ss));
        denseRatio(ss) = max(Yfit)/min(Yfit);
        if denseRatio(ss)>30
            loglog(supportDeg,Yfit,'-b');
        else
            loglog(supportDeg,Yfit,'-k');
        end
    end
    ylabel('log_1_0 density [cones/deg^2]')
    xlabel('Eccentricity [deg]');
    title([meridianLabels{meridianIdx} ' meridian subject fits']);
end
plotFileName = fullfile(sourceDir,'figures','Fig08_individualSubjectNasalMerdianFits.pdf');
saveas(gcf,plotFileName);


%% Plot individual subject fits
for ss=1:length(subNames)
    
    if missingSplit(ii)
        continue
    end
    
    Y = squeeze(dataMat(:,:,ss));
    Yfit = squeeze(YfitSet(:,:,ss));
    p = pSet(:,ss);

    figHandle = figure();
    
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
    
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_polarModelFit.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);

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

