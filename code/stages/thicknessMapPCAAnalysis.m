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
p.addParameter('nCoeff',5,@isscalar);
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
%         fileName = fullfile(rawDirList(ss).folder,rawDirList(ss).name,[rawDirList(ss).name '_' p.Results.layerSetLabels{ii} '_volumeMap.mat']);
%         load(fileName,'volumeMap_mmCubedDegSquared');
%         thisMap = volumeMap_mmCubedDegSquared;
       fileName = fullfile(rawSubjectList(ss).folder,rawSubjectList(ss).name);
       load(fileName,'averageMaps');
       thisMap = averageMaps.(p.Results.layerSetLabels{ii});
       % Convert to mm and then to cubic mm of tissue per degree
       thisMap = thisMap./1000 .* mmPerDeg(axialLength(ss))^2;
        if ss==1
            imageSize = size(thisMap);
            avgMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
        end
        % Convert the map from thickness to volume per degree^2
        avgMapBySubject(ss,:,:)=thisMap;

        % Add this ID to the table
        pcaResultsTable.AOSO_ID(ss) = str2num(rawDirList(ss).name);
    end
        
    % Determine which points are not nans in any subject
    observedIdx = ~isnan(squeeze(sum(avgMapBySubject,1)));
    
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
    
    % Find the two subjects with the values close to the median for the
    % first coefficient, but maximally different from each other on the
    % second coefficient
    medianThickness = median(score(:,1));
    [~,idxThickOrder] = sort(abs(score(:,1)-medianThickness),'ascend');
    rankThickOrder(idxThickOrder)=1:length(idxThickOrder);

    [~,idxLowCoeff] = sort(score(:,2),'ascend');
    rankLowCoeff2(idxLowCoeff)=1:length(idxLowCoeff);
    [~,idxHighCoeff] = sort(score(:,2),'descend');
    rankHighCoeff2(idxHighCoeff)=1:length(idxHighCoeff);
    [~,idxLowCoeff] = sort(score(:,3),'ascend');
    rankLowCoeff3(idxLowCoeff)=1:length(idxLowCoeff);
    [~,idxHighCoeff] = sort(score(:,3),'descend');
    rankHighCoeff3(idxHighCoeff)=1:length(idxHighCoeff);
    
    rankThickLow = mean([rankThickOrder; rankLowCoeff2; rankLowCoeff3]);
    rankThickHigh = mean([rankThickOrder; rankHighCoeff2; rankHighCoeff3]);
    [~,subjectIdx(1)] = min(rankThickLow);
    [~,subjectIdx(2)] = min(rankThickHigh);
    
    
    % Show a plot of the explanatory power of the components
    if p.Results.showPlots
        figure
        plot(explained)
        hold on
        plot(explained,'or')
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'varianceExplained.pdf');
            print(gcf,filename,'-dpdf');
        end
    end
    
    % Create a 3D plot of the scores
    if p.Results.showPlots
        figure
        plot3(score(:,1),score(:,2),score(:,3),'*k')
        hold on
        plot3(score(subjectIdx(1),1),score(subjectIdx(1),2),score(subjectIdx(1),3),'*r')
        plot3(score(subjectIdx(2),1),score(subjectIdx(2),2),score(subjectIdx(2),3),'*b')
        axis equal
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
            plot(axialLength,score(:,jj),'ok')
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
    
    % Now show images of the first, second, and third PCA components
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
        figA = figure();
        figB = figure();
        zMax = 0.01;
            supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1));
            [xMesh,yMesh] = meshgrid(supportDeg,supportDeg);

        for jj = 1:length(subjectIdx)
            subjectID = rawDirList(jj).name;
            
            figure(figA);
            subplot(length(subjectIdx),3,1+3*(jj-1))
            theMap = squeeze(avgMapBySubject(subjectIdx(jj),:,:));
            mesh(xMesh,yMesh,theMap);
            caxis([0 zMax]);
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            zlim([0 zMax]);
            title(['Subject ' subjectID] )
            axis square
            
            subplot(length(subjectIdx),3,2+3*(jj-1))
            reconMap = nan(imageSize(1),imageSize(2));
            reconMap(:)=XfitFit(subjectIdx(jj),:);
            mesh(xMesh,yMesh,reconMap);
            caxis([0 zMax]);
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            zlim([0 zMax]);
            title('PCA reconstructed')
            axis square

            subplot(length(subjectIdx),3,3+3*(jj-1))
            mesh(xMesh,yMesh,reconMap-theMap);
            caxis([-zMax/2 zMax/2]);
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            zlim([-zMax/2 zMax/2]);
            title('Error')
            axis square
            
            figure(figB);
            plot(supportDeg,reconMap(imageSize(1)/2,:));
            hold on
            xlim([-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent]);
            ylim([0 zMax]);
            title('Horizontal meridian')
            axis square
        end
        if ~isempty(p.Results.figSaveDir)
            filename = fullfile(p.Results.figSaveDir,'pcaDataReconMap');
            vecrast(figA, filename, 600, 'top', 'pdf')
            filename = fullfile(p.Results.figSaveDir,'pcaDataReconMeridian.pdf');
            print(figB,filename,'-dpdf');
        end

    end
    
    % Examine the correlation of PCA scores with anatomical measures
    
    % Join the anat measures to the subject data table, and then to the PCA
    % scores pcaResultsTable
    comboTable = join(pcaResultsTable,subjectTable,'Keys','AOSO_ID');
    comboTable = join(anatMeasuresTable,comboTable,'Keys','AOSO_ID');

    chiasmRelative = comboTable.Optic_Chiasm./comboTable.SupraTentorialVol;
    v1AreaRelative = (comboTable.lh_lh_v1_noah_template_label_area+comboTable.rh_rh_v1_noah_template_label_area)./(comboTable.lh_WhiteSurfArea_area+comboTable.rh_WhiteSurfArea_area);
    v1ThickRelative = (comboTable.lh_lh_v1_noah_template_label_thickness+comboTable.rh_rh_v1_noah_template_label_thickness)./(comboTable.lh_lh_cortex_label_thickness+comboTable.rh_rh_cortex_label_thickness);

    X = [comboTable.PCA1 comboTable.PCA2 -comboTable.PCA2.^2 comboTable.PCA3 -comboTable.PCA3.^2 comboTable.Axial_Length_average];
    X = X-mean(X);
    
    lm = fitlm(X,comboTable.Optic_Chiasm)    
    lm = fitlm(X,comboTable.LGNJacobian)    
    lm = fitlm(X,v1AreaRelative)
    lm = fitlm(X,v1ThickRelative)

    fopo=1;
end

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