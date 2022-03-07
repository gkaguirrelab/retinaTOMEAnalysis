%% polarModelIllustrations
% This routine creates plots that describe the elements of the polar model
% of cone density

clear

% Set up some constants
supportLengthDeg = 1799;
imRdim = (supportLengthDeg+1)/2;
supportDegDelta = 0.00773;
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLengthDeg-1);
maxSupportDeg = 15;
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];
supportPA = linspace(0,360,supportLengthDeg);
polarRatio = (supportLengthDeg+1)/360;

% Load the individual subject fit
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsDeg.mat');
load(individualFitFile,'pSet','YfitSet','fValSet','RSquaredSet','polarMultiplierSet','dataMatDeg','subNames')

% Create the X and P support, the mean data, and the fit
X = repmat(supportDeg,supportLengthDeg,1);
P = repmat(linspace(0,360,supportLengthDeg)',1,supportLengthDeg);
Y = nanmean(dataMatDeg,3);
p = mean(pSet,2);
w = sum(~isnan(dataMatDeg),3);

% Create the Yfit
Yfit = coneDensityModel(X,P,maxSupportDeg,p);



% Illustrate the model surface components
figure
for cc = 1:4
    subplot(2,2,cc)
    switch cc
        case 1
            Ymodel = cosd(P+p(5));
        case 2
            Ymodel = sind(P+p(5));
        case 3
            Ymodel = cosd(P.*2+p(13));
        case 4
            Ymodel = cosd(P.*4+p(17));
    end
    surf(X,P,Ymodel,'FaceAlpha',0.5,'EdgeColor','none');
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [deg]');
    zlabel('Density modulation');
    view(45,15)
end
plotFileName = fullfile(sourceDir,'figures','Fig0X_polarModulationModelDeg.pdf');
saveas(gcf,plotFileName);

% Illustrate the model surface components
figure
for cc = 1:4
    subplot(2,2,cc)
    switch cc
        case 1
            g = p(6).*gampdf(X,p(7),p(8))./max(gampdf(0:0.01:maxSupportDeg,p(7),p(8)));
            Ymodel = g.*cosd(P+p(5)+180);
        case 2
            g = p(10).*gampdf(X,p(11),p(12))./max(gampdf(0:0.01:maxSupportDeg,p(11),p(12)));
            Ymodel = g.*sind(P+p(9)+180);
        case 3
            g = p(14).*gampdf(X,p(15),p(16))./max(gampdf(0:0.01:maxSupportDeg,p(15),p(16)));
            Ymodel = g.*cosd(P.*2+p(13));
        case 4
            g = p(18).*gampdf(X,p(19),p(20))./max(gampdf(0:0.01:maxSupportDeg,p(19),p(20)));
            Ymodel = g.*cosd(P.*4+p(17));
        case 5
            g = p(22).*gampdf(X,p(23),p(24))./max(gampdf(0:0.01:maxSupportDeg,p(23),p(24)));
            Ymodel = g.*sind(P+p(21)+180);
    end
    surf(X,P,Ymodel,'FaceAlpha',0.5,'EdgeColor','none');
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [deg]');
    zlabel('Density modulation');
    zlim([-0.1 0.1]);
    view(45,15)
end
plotFileName = fullfile(sourceDir,'figures','Fig0X_polarModulationModelWithGammaDeg.pdf');
saveas(gcf,plotFileName);

% Illustrate the model surface components
figure
g = p(6).*gampdf(X,p(7),p(8))./max(gampdf(0:0.01:maxSupportDeg,p(7),p(8)));
Ymodel = g.*cosd(P+p(5)+180);
g = p(10).*gampdf(X,p(11),p(12))./max(gampdf(0:0.01:maxSupportDeg,p(11),p(12)));
Ymodel = Ymodel+g.*sind(P+p(9)+180);
g = p(14).*gampdf(X,p(15),p(16))./max(gampdf(0:0.01:maxSupportDeg,p(15),p(16)));
Ymodel = Ymodel+g.*cosd(P.*2+p(13));
g = p(18).*gampdf(X,p(19),p(20))./max(gampdf(0:0.01:maxSupportDeg,p(19),p(20)));
Ymodel = Ymodel+g.*cosd(P.*4+p(17));
supportPA = linspace(0,360,supportLengthDeg);
contourf(supportDeg,supportPA,Ymodel,15,'LineWidth',2)
map = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';
colormap(map)

%imagesc(supportDeg,supportPA,Ymodel)
%colorbar
ax = gca;
ax.YDir = 'normal';
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density modulation');
plotFileName = fullfile(sourceDir,'figures','Fig0X_entirePolarModulationDeg.pdf');
saveas(gcf,plotFileName);




figure
legendLabels = {'cos1','sin1','cos2','cos4'};
nFourier = length(legendLabels);
for gg = 1:nFourier
    ph = p((gg-1)*4+5);
    f1 = p((gg-1)*4+6);
    f2 = p((gg-1)*4+7);
    f3 = p((gg-1)*4+8);
    g = f1.*gampdf(supportDeg,f2,f3)./max(gampdf(0:0.01:maxSupportDeg,f2,f3));
    plot(supportDeg,g)
    hold on
end
legend(legendLabels(1:nFourier));
xlabel('Eccentricity [deg]');
ylabel('Modulation [proportion]');
plotFileName = fullfile(sourceDir,'figures','Fig0X_gammaPDFWeightFunctionsDeg.pdf');
saveas(gcf,plotFileName);

