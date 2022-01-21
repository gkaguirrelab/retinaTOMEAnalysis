clear

% Set up some constants
supportLengthMm = 2001;
imRdim = (supportLengthMm+1)/2;
supportMmDelta = 0.0025;
supportMm = 0:supportMmDelta:supportMmDelta*(supportLengthMm-1);
maxSupportMm = 5;
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];
supportPA = linspace(0,360,supportLengthMm);
polarRatio = (supportLengthMm+1)/360;

% Load the individual subject fit
sourceDir = '/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images/densityAnalysis/';
individualFitFile = fullfile(sourceDir,'individualSubjectFitsMm.mat');
load(individualFitFile,'pSet','YfitSet','fValSet','RSquaredSet','polarMultiplierSet','dataMatMm','subNames')

% How many subjects?
nSubs = length(subNames);


%% Plot individual subject diagnostic plots

for ss=1:length(subNames)
    
    Y = squeeze(dataMatMm(:,:,ss));
    Yfit = squeeze(YfitSet(:,:,ss));
    p = pSet(:,ss);
    
    % Mean polar data and model fit
    figHandle = figure();
    X = repmat(supportMm,supportLengthMm,1);
    P = repmat(linspace(0,360,supportLengthMm)',1,supportLengthMm);
    surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none');
    hold on
    plot3(X(:),P(:),Y(:),'.k')
    yticks(meridianAngles);
    yticklabels(meridianLabels);
    xlabel('Eccentricity [mm]');
    zlabel('Density [cones/deg^2]');
    view(45,15)
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_polarSurfaceMm.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);
        
    % Polar bands and model fit
    figHandle = figure();    
    for ii = [0.03125   0.0625    0.1250    0.2500    0.5000    1.0000    2.0000]
        idx = find(supportMm>=ii,1);
        if isempty(idx)
            density = coneDensityModel(ii,0,maxsupportMm,pSet(:,ss));
            semilogy(repmat(density,1,size(dataMatMm,1)),'-r');
            text(300*polarRatio,density,sprintf('%2.1f mm',ii));
        else
            semilogy(Y(:,idx),'.');
            hold on
            semilogy(Yfit(:,idx),'-r');
            text(300*polarRatio,Yfit(round(45*polarRatio),idx),sprintf('%2.1f°',supportMm(idx)));
        end
    end
    xticks(meridianAngles*polarRatio);
    xticklabels(meridianLabels);
    ylim([10^2,2*10^4]);
    ylabel('log_1_0 density [cones/deg^2]')
    
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_polarModelFitMm.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);
    
    % Plot by meridian
    figHandle = figure();
    
    for mm=1:4
        subplot(2,2,mm)
        plot(supportMm,Y(round((meridianAngles(mm))*polarRatio+1),:),'.k');
        hold on
        plot(supportMm,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
        xlabel('Eccentricity [mm]');
        ylabel('Density [cones/deg^2]');
        title(meridianLabels{mm});
    end
    
    plotFileName = fullfile(sourceDir,'figures','subjects',[subNames{ss} '_MeridianModelFitMm.pdf']);
    saveas(figHandle,plotFileName);
    close(figHandle);    
    
end