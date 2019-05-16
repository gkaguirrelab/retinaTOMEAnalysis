function thicknessMapPCAAnalysis(dataDir, varargin)
% Do some analysis
%
% Description:
%   Foo
%


%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('dataDir',@ischar);

% Optional analysis params
p.addParameter('layerSetLabels',{'RGCIPL'},@iscell);
p.addParameter('showPlots',true,@islogical);
p.addParameter('figSaveDir','/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Aguirre_2019_rgcCorticalAnatomy/VSS2019/raw figures/pca',@ischar);
p.addParameter('makeFitMaps',false,@islogical);
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);
p.addParameter('anatMeasuresFileName',fullfile(getpref('retinaTOMEAnalysis','projectBaseDir'),'data','visualPathwayAnatMeasures.xlsx'),@ischar);
p.addParameter('nCoeff',4,@isscalar);
p.addParameter('octRadialDegreesVisualExtent',15,@isscalar);


%% Parse and check the parameters
p.parse(dataDir, varargin{:});

% Obtain a list of subjects
rawSubjectList = dir(fullfile(dataDir,'*/*.mat'));
rawDirList = dir(fullfile(dataDir,'1*'));
nSubs = length(rawDirList);

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);
axialLength = subjectTable.Axial_Length_average;

% Load the anat measures
opts = detectImportOptions(p.Results.anatMeasuresFileName);
anatMeasuresTable = readtable(p.Results.anatMeasuresFileName, opts);

% The number of PCA coefficients we will retain
nCoeff = p.Results.nCoeff;
variableNames = {'AOSO_ID'};
for kk = 1:nCoeff
    variableNames = [variableNames {['PCA' num2str(kk)]}];
end

% This is the mmPerDeg at the ellipsoidal pole of the vitreous chamber
mmPerDeg = @(axialLength) (0.0165.*axialLength)-0.1070;

% Loop over layer sets
for ii = 1:length(p.Results.layerSetLabels)
    
    % Create a table to hold the PCA results
    pcaResultsTable = table('Size',[nSubs,nCoeff+1],'VariableTypes',repmat({'double'},1,nCoeff+1),'VariableNames',variableNames);
    
    % Loop over subjects and load the maps
    for ss = 1:nSubs
        fileName = fullfile(rawSubjectList(ss).folder,rawSubjectList(ss).name);
        load(fileName,'averageMaps');
        thisMap = averageMaps.(p.Results.layerSetLabels{ii});
        % Convert from microns to mm
        thisMap = thisMap./1000;
        % If this is the first subject, set up some variables
        if ss==1
            imageSize = size(thisMap);
            avgMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
            rawMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
        end
        % Store the raw map
        rawMapBySubject(ss,:,:)=thisMap;
        % Convert to cubic mm of tissue per degree
        thisMap = thisMap .* mmPerDeg(axialLength(ss))^2;
        % Store the map
        avgMapBySubject(ss,:,:)=thisMap;
        % Add this ID to the table
        pcaResultsTable.AOSO_ID(ss) = str2double(rawDirList(ss).name);
    end
    
    % Determine which points are not nans in any subject
    observedIdx = ~isnan(squeeze(sum(avgMapBySubject,1)));
    
    % Determine the relationship between mean thickness and axial length
    rawVecBySubject = reshape(rawMapBySubject,50,prod(imageSize));
    meanThickBySubject = mean(rawVecBySubject(:,observedIdx(:)),2);
    if p.Results.showPlots
        figure
        plot(axialLength,meanThickBySubject,'.k');
        hold on
        pFit = polyfit(axialLength,meanThickBySubject,1);
        yfit = pFit(1)*axialLength+pFit(2);
        plot(axialLength,yfit,'r-');
        xlabel('axial length [mm]');
        xlim([20 28]);
        ylabel('mean thickness [microns]');
        ylim([0 0.075]);
        rVal = corr2(axialLength,meanThickBySubject);
        title(['Pearson r = ' num2str(rVal)]);
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'meanThickByAxialLength.pdf');
            print(gcf,filename,'-dpdf');
        end
    end
    
    
    % Conduct a PCA analysis across subjects, averaged over eyes. First,
    % Convert the non-nan portion of the image to a vector.
    X = [];
    for ss=1:nSubs
        tmp = squeeze(avgMapBySubject(ss,:,:));
        X(ss,:) = tmp(observedIdx);
    end
    
    % Calc the PCA
    [coeff,score,~,~,explained,mu] = pca(X,'Centered',false);
    
    % To simplify subsequent interpretation, reverse the sign of the
    % second PCA component
    score(:,2) = -score(:,2);
    coeff(:,2) = -coeff(:,2);
    
    % Store the results in the data table
    for kk=1:nCoeff
        pcaResultsTable.(['PCA' num2str(kk)]) = score(:,kk);
    end
    
    % Join the anat measures to the subject data table, and then to the PCA
    % scores pcaResultsTable
    comboTable = join(pcaResultsTable,subjectTable,'Keys','AOSO_ID');
    comboTable = join(anatMeasuresTable,comboTable,'Keys','AOSO_ID');
    
    % Derive some anatomical measures
    brainMeasureNames = {'chiasmAbsolute','lgnRelative','v1AreaRelative','v1AreaAbsolute','v1ThickRelative','v1ThickAbsolute'};
    chiasmAbsolute = comboTable.Optic_Chiasm;
    lgnRelative = comboTable.LGNJacobian;
    v1AreaRelative = 100.*((comboTable.lh_lh_v1_noah_template_label_area+comboTable.rh_rh_v1_noah_template_label_area)./(comboTable.lh_WhiteSurfArea_area+comboTable.rh_WhiteSurfArea_area));
    v1AreaAbsolute = comboTable.lh_lh_v1_noah_template_label_area+comboTable.rh_rh_v1_noah_template_label_area;
    v1ThickRelative = 100.*((comboTable.lh_lh_v1_noah_template_label_thickness+comboTable.rh_rh_v1_noah_template_label_thickness)./(comboTable.lh_lh_cortex_label_thickness+comboTable.rh_rh_cortex_label_thickness));
    v1ThickAbsolute = (comboTable.lh_lh_v1_noah_template_label_thickness+comboTable.rh_rh_v1_noah_template_label_thickness)./2;
    
    % Create a design matrix composed of the PCA scores.
    X = [comboTable.PCA1 comboTable.PCA2 comboTable.PCA3 comboTable.Axial_Length_average];
    X = X-mean(X);
    
    % Examine the relation between PCA scores and brain measures.
    plotChoices = [1 3];
    plotChoicesYLabels = {'optic chiasm volume [mm^3]','V1 relative area [%]'};
    plotChoicesYLims = {[100 500],[2 4]};
    plotChoicesZMax = {0.015,0.010};
    plotChoicesSymbols = {'x','o'};
    subjectSets = zeros(length(plotChoices),2);
    fprintf('\nRobust linear regression, OCT PCA 3 components and axial length, upon brain measures:\n\n');
    for kk = 1:length(brainMeasureNames)
        thisVarName = brainMeasureNames{kk};
        thisVar = eval(thisVarName);
        outline=['   ' thisVarName ' -  '];
        lm = fitlm(X,thisVar,'RobustOpts','on');
        outline = [outline 'p-value: ' num2str(lm.coefTest) '  '];
        outline = [outline 'adjusted R^2: ' num2str(lm.Rsquared.Adjusted) '  '];
        outline = [outline 'coeff p-vals: ' num2str(lm.Coefficients.pValue(2:end)') '  '];
        fprintf([outline '\n']);
        if any(plotChoices==kk)
            subjectSets(1,:) = [14 2];
            subjectSets(2,:) = [36 7];
            % Select two subjects who are the extreme ends of the fitted
            % values for this relationship
%            [~,subjectSets(plotChoices==kk,1)]=min(lm.Fitted);
%            [~,subjectSets(plotChoices==kk,2)]=max(lm.Fitted);
            % Plot the linear model plot
            if p.Results.showPlots
                figure
                lm.plot
                hold on
%                plot(lm.Fitted(subjectSets(plotChoices==kk,1)),thisVar(subjectSets(plotChoices==kk,1)),'or')
%                plot(lm.Fitted(subjectSets(plotChoices==kk,2)),thisVar(subjectSets(plotChoices==kk,2)),'ob')
                ylabel(plotChoicesYLabels{plotChoices==kk});
                ylim(plotChoicesYLims{plotChoices==kk});
                if ~isempty(p.Results.figSaveDir)
                    filename = fullfile(p.Results.figSaveDir,['robustRegression_' thisVarName '.pdf']);
                    print(gcf,filename,'-dpdf');
                end
            end
        end
    end
    fprintf('\n');
    
    % Report variance explained
    explained = [nan; explained(2:end)./sum(explained(2:end))];
    outline=sprintf('Variance explained by first %d coefficients: %2.2f \n',nCoeff,nansum(explained(1:nCoeff)));
    fprintf(outline);
    
    % Report the mean values of the fitted coefficients
    outline=sprintf('Mean thickness of first four coefficients for this population: [%2.4f, %2.4f, %2.4f, %2.4f]\n',...
        mean(mean(score(:,1)*coeff(:,1)')),...
        mean(mean(score(:,2)*coeff(:,2)')),...
        mean(mean(score(:,3)*coeff(:,3)')),...
        mean(mean(score(:,4)*coeff(:,4)')));
    fprintf(outline);
    
    % Obtain the fit to the data
    Xfit = score(:,1:nCoeff)*coeff(:,1:nCoeff)'+(repmat(mu,nSubs,1));
           
    % Show a plot of the explanatory power of the components
    if p.Results.showPlots
        figure
        plot(explained)
        hold on
        plot(explained,'.r')
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'varianceExplained.pdf');
            print(gcf,filename,'-dpdf');
        end
    end
    
    % Create a 3D plot of the scores
    if p.Results.showPlots
        figure
        plot3(score(:,1),score(:,2),score(:,3),'.k')
        hold on
        for xx=1:size(subjectSets,1)
            plot3(score(subjectSets(xx,:),1),score(subjectSets(xx,:),2),score(subjectSets(xx,:),3),'-k')
            plot3(score(subjectSets(xx,1),1),score(subjectSets(xx,1),2),score(subjectSets(xx,1),3),[plotChoicesSymbols{xx} 'r'])
            plot3(score(subjectSets(xx,2),1),score(subjectSets(xx,2),2),score(subjectSets(xx,2),3),[plotChoicesSymbols{xx} 'b'])
        end
        axis square
        xlabel('PCA1');
        ylabel('PCA2');
        zlabel('PCA3');
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'subjectScores.pdf');
            print(gcf,filename,'-dpdf');
        end
    end
    
    % Plot the relationship between the scores and axial length
    if p.Results.showPlots
        figure
        for jj=1:nCoeff
            subplot(2,2,jj);
            plot(axialLength,score(:,jj),'.k')
            corrVal = corr2(squeeze(score(:,jj)),axialLength);
            title(['R = ' num2str(corrVal)]);
        end
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'scoresVsAxialLength.pdf');
            print(gcf,filename,'-dpdf');
        end
    end
    
    % Create interpolated versions of the PCA coefficient maps
    filename = fullfile(p.Results.figSaveDir,'pcaCoeffMapsFit.mat');
    if exist(filename,'file')
        load(filename,'fitCoeffMaps','fitCoeff')
    else
        for jj = 1:nCoeff
            coeffMap = nan(imageSize(1),imageSize(2));
            coeffMap(observedIdx)=coeff(:,jj);
            fitMap = thinSplineFitMap(coeffMap, p.Results.octRadialDegreesVisualExtent);
            fitCoeff(:,jj) = fitMap(:);
            fitCoeffMaps(jj,:,:)=fitMap;
        end
        save(filename,'fitCoeffMaps','fitCoeff')
    end
    XfitFit = score(:,1:nCoeff)*fitCoeff(:,1:nCoeff)';
    if p.Results.showPlots
        figure
        for jj=1:nCoeff
            subplot(2,2,jj);
            supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1));
            [xMesh,yMesh] = meshgrid(supportDeg,supportDeg);
            mesh(xMesh,yMesh,squeeze(fitCoeffMaps(jj,:,:)));
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            axis square
            title('PCA coefficients (fit)')
            set(gca,'color','none')
        end
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'pcaCoeffMapsFit');
            vecrast(gcf, filename, 600, 'top', 'pdf')
        end
    end
    
    % Now show images of the PCA components
    coeffMap = nan(imageSize(1),imageSize(2));    
    if p.Results.showPlots
        figure
        for jj=1:nCoeff
            subplot(2,2,jj);
            supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1));
            [xMesh,yMesh] = meshgrid(supportDeg,supportDeg);
            coeffMap(observedIdx)=coeff(:,jj);
            mesh(xMesh,yMesh,coeffMap);
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            axis square
            title('PCA coefficients (raw)')
            set(gca,'color','none')
        end
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'pcaCoeffMapsRaw');
            vecrast(gcf, filename, 600, 'top', 'pdf')
        end
    end
    
    % Show some PCA reconstructed OCT scans vs. the actual data
    if p.Results.showPlots
        supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1));
        [xMesh,yMesh] = meshgrid(supportDeg,supportDeg);
        
        for xx = 1:size(subjectSets,1)
            
            figA = figure();
            figB = figure();

            zMax = plotChoicesZMax{xx};
            
            for jj = 1:size(subjectSets,2)
                subjectID = rawDirList(subjectSets(xx,jj)).name;
                
                figure(figA);
                subplot(size(subjectSets,2),3,1+3*(jj-1))
                theMap = squeeze(avgMapBySubject(subjectSets(xx,jj),:,:));
                mesh(xMesh,yMesh,theMap);
                caxis([0 zMax]);
                xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                zlim([0 zMax]);
                zlabel('tissue volume [mm3/deg2]');
                title(['Subject ' subjectID] )
                axis square
                
                subplot(size(subjectSets,2),3,2+3*(jj-1))
                reconMap = nan(imageSize(1),imageSize(2));
                reconMap(:)=XfitFit(subjectSets(xx,jj),:);
                mesh(xMesh,yMesh,reconMap);
                caxis([0 zMax]);
                xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                zlim([0 zMax]);
                zlabel('tissue volume [mm3/deg2]');
                title('PCA reconstructed')
                axis square
                
                subplot(size(subjectSets,2),3,3+3*(jj-1))
                mesh(xMesh,yMesh,reconMap-theMap);
                caxis([-zMax/2 zMax/2]);
                xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                zlabel('tissue volume [mm3/deg2]');
                zlim([-zMax/2 zMax/2]);
                title('Error')
                axis square
                
                figure(figB);
                plot(supportDeg,reconMap(imageSize(1)/2,:));
                hold on
                xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
                ylim([0 zMax]);
                zlabel('tissue volume [mm3/deg2]');
                title('Horizontal meridian')
                axis square
            end
            if ~isempty(p.Results.figSaveDir)
                filename = fullfile(p.Results.figSaveDir,['pcaDataReconMap_MaxMinSubs_' brainMeasureNames{plotChoices(xx)}]);
                vecrast(figA, filename, 600, 'top', 'pdf')
                filename = fullfile(p.Results.figSaveDir,['pcaDataReconMeridian_' brainMeasureNames{plotChoices(xx)} '.pdf']);
                print(figB,filename,'-dpdf');
                close(figA);
                close(figB);
            end
        end
    end
    
end % Loop over layer sets

end % Main function


function fitMap = thinSplineFitMap(rawMap, radialDegrees)

% Downsample the rawMap by a factor of 12
scaleDown = 12;
imageSize = size(rawMap);
supportDeg = linspace(-radialDegrees,radialDegrees,imageSize(1)/scaleDown);
downSampMap = imresize(rawMap,1/scaleDown);

% Silence a warning that occurs regarding nans in the maps
warningState = warning;
warning('off','curvefit:prepareFittingData:removingNaNAndInf');

% Prepare the surface and perform the fit
[xo,yo,zo]=prepareSurfaceData(supportDeg,supportDeg,downSampMap);
sf = fit([xo, yo],zo,'thinplateinterp');

% Restore the warning state
warning(warningState);

% Create the interpolated map just within the radius
supportDeg = linspace(-radialDegrees,radialDegrees,imageSize(1));
[xi,yi]=meshgrid(supportDeg,supportDeg);
inRangeIdx = sqrt(xi.^2+yi.^2)<=radialDegrees;
fitCoeffVals = sf(xi(inRangeIdx),yi(inRangeIdx));
fitMap = nan(imageSize(1),imageSize(2));
fitMap(inRangeIdx)=fitCoeffVals;

end