%% groupModelFitFigures


% Set up some constants
supportLength = 1799;
imRdim = (supportLength+1)/2;
supportDegDelta = 0.00773;
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLength-1);
maxSupportDeg = 15;
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];
supportPA = linspace(0,360,supportLength);

% Load the individual subject fit
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFits.mat');
load(individualFitFile,'pSet','YfitSet','fValSet','RSquaredSet','polarMultiplierSet','dataMat','subNames')

% Create the X and P support, the mean data, and the fit
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);
Y = nanmean(dataMat,3);
p = mean(pSet,2);
w = sum(~isnan(dataMat),3);

% Create the Yfit
Yfit = coneDensityModel(X,P,maxSupportDeg,p);

% Extract a couple of costants
nSubs = size(pSet,2);
polarRatio = (size(dataMat,1)+1)/360;

% Some machinery for bootstrap estimation of fit SEM
%{
nSubs = size(pSet,2);
nBoots = 10;
bootYfit = nan(supportLength,supportLength,nBoots);
for bb = 1:nBoots
    boot = randsample(nSubs,nSubs,true);
    p = mean(pSet(:,boot),2);
    bootYfit(:,:,bb) = coneDensityModel(X,P,maxSupportDeg,p);
end
semYfit = std(bootYfit,0,3);
Yfit = mean(bootYfit,3);
%}


%% Figures describing the mean model fit

% Plot the mean polar data and model fit
figHandle = figure;

contourf(X,P,Yfit,logspace(log10(500),log10(15000),15),'LineWidth',2)
map = [ logspace(log10(0.5),log10(1),255); logspace(log10(0.5),log10(0.1),255); logspace(log10(0.5),log10(0.1),255)]';
colormap(map)
% 
% 
% sHandle = surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none');
% hold on
% %plot3(X(:),P(:),Y(:),'.k')
% yticks(meridianAngles);
% yticklabels(meridianLabels);
% xlabel('Eccentricity [deg]');
% zlabel('Density [cones/deg^2]');
% % Add topographic lines
% ridgeHandles = [];
% for targetDensity = 500:500:15000
%     eccenIdx = nan(supportLength,1);
%     for ii=1:supportLength
%         [~,eccenIdx(ii)]=min(abs(Yfit(ii,:)-targetDensity));
%     end
%     ridgeHandles(end+1) = plot3(supportDeg(eccenIdx),supportPA,repmat(targetDensity,1,supportLength),'-k','LineWidth',1.5);
% end
% view(10,10)
% lighting gouraud
% lightangle(gca,10,10)
% lightangle(gca,10,10)
% lightangle(gca,10,10)
% lightangle(gca,0,0)
% lightangle(gca,0,0)
% lightangle(gca,0,0)
% plotFileName = fullfile(sourceDir,'figures','Fig01_meanModelFitPolar');
% alignHandle = plot3(15,360,16000,'xm');
% zlim([0 16000]);
% 
% % Save figure with vector components
% hidem(sHandle);
% set(figHandle,'color','white');
% fileName = [plotFileName '.pdf'];
% export_fig(figHandle,fileName,'-Painters');
% 
% % Save figure with rendered components
% axis off
% showm(sHandle);
% hidem(ridgeHandles);
% set(figHandle,'color','none');
% fileName = [plotFileName '.png'];
% export_fig(figHandle,fileName,'-r1200','-opengl');
% 




% Mean Cartesian data and model fit
figure
cartSupportDeg=linspace(-max(supportDeg),max(supportDeg),imRdim);
[cartXDeg,cartYDeg] = meshgrid(cartSupportDeg,cartSupportDeg);
cartYfit = convertPolarMapToImageMap(Yfit,'imRdim',imRdim);
cartYfit(cartYfit<min([Y(:); Yfit(:)]))=nan;
cartY = convertPolarMapToImageMap(Y,'imRdim',imRdim);
cartY(cartY<min([Y(:); Yfit(:)]))=nan;



contourf(cartSupportDeg,cartSupportDeg,cartYfit,logspace(log10(500),log10(15000),15),'LineWidth',2)
map = [ logspace(log10(0.5),log10(1),255); logspace(log10(0.5),log10(0.1),255); logspace(log10(0.5),log10(0.1),255)]';
colormap(map)

xlabel('Eccentricity [deg]');
ylabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');

axis off

plotFileName = fullfile(sourceDir,'figures','Fig02_meanModelFitCartesian.pdf');
saveas(gcf,plotFileName);

% Polar angle density variation
figure
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
plotFileName = fullfile(sourceDir,'figures','Fig03_densityVariationByPolarAngle.pdf');
saveas(gcf,plotFileName);

% Density and model by meridian
figure
for mm=1:4
    subplot(2,2,mm)
    yData = Y(round((meridianAngles(mm))*polarRatio+1),:);
    % Mask the fovea data points. When we fit the group model, we don't
    % include these points as there is a bias towards lower vaues.
    yData(supportDeg<0.3) = nan;
    plot(supportDeg,yData,'.k');
    hold on
    plot(supportDeg,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
    xlabel('Eccentricity [deg]');
    ylabel('Density [cones/deg^2]');
    title(meridianLabels{mm});
end
plotFileName = fullfile(sourceDir,'figures','Fig04_meanModelFitByMeridian.pdf');
saveas(gcf,plotFileName);

% Model fit by meridian
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

% Individual density functions by meridian
figure
for meridianIdx = 1:4
    subplot(2,2,meridianIdx)

    % Add the Curcio profile along the specified meridian
    CurcioFitConeDensitySqDegVisual = getSplineFitToConeDensitySqDegVisual(meridianAngles(meridianIdx));
    loglog(supportDeg,CurcioFitConeDensitySqDegVisual(supportDeg),'-r','LineWidth',5);

    % Add the individual subject fits
    hold on
    for ss = 1:nSubs
        Yfit = squeeze(YfitSet(round((meridianAngles(meridianIdx))*polarRatio+1),:,ss));
        loglog(supportDeg,Yfit,'-k');
    end
    ylabel('density [cones/deg^2]')
    xlabel('eccentricity [deg]');
    title([meridianLabels{meridianIdx} ' meridian subject fits']);
end
plotFileName = fullfile(sourceDir,'figures','Fig06_individualDensityFcnsByMeridian.pdf');
saveas(gcf,plotFileName);

% Loop over subjects and calculate the integrated cone count
for ss = 1:size(dataMat,3)    
    myFunc = @(ec,pa) coneDensityModel(ec,pa,maxSupportDeg,pSet(:,ss));
    count(ss) = integral2(myFunc,0,10,0,360);
end

% Map of weights (number of subjects) at each polar map location
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



figure
histogram(count)
xlabel('integrated cone density')
ylabel('counts [subjects]')


figure
rSquareTitles = {'RSquaredFull','RSquaredExponentialOnly','RSquaredPolarOnly','RSquaredPolarResiduals'};
for ii=1:4
    subplot(2,2,ii)
    histogram(RSquaredSet(ii,:))
xlabel('R squared')
ylabel('counts [subjects]')
title(rSquareTitles{ii})
end

figure
histogram(polarMultiplierSet)
xlabel('polar model multiplier')
ylabel('counts [subjects]')
