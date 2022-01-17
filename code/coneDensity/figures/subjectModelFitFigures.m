

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