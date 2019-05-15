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
p.addParameter('subjectTableFileName',fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'),'TOME_subject','TOME-AOSO_SubjectInfo.xlsx'),@ischar);
p.addParameter('anatMeasuresFileName',fullfile(getpref('retinaTOMEAnalysis','projectBaseDir'),'data','visualPathwayAnatMeasures.xlsx'),@ischar);
p.addParameter('nCoeff',4,@isscalar);
p.addParameter('octRadialDegreesVisualExtent',15,@isscalar);


%% Parse and check the parameters
p.parse(dataDir, varargin{:});

% Obtain a list of subjects
rawSubjectList = dir(fullfile(dataDir,'*/*.mat'));
rawDirList = dir(fullfile(dataDir,'1*'));
nSubs = length(rawSubjectList);

% Load the subject data table
opts = detectImportOptions(p.Results.subjectTableFileName);
subjectTable = readtable(p.Results.subjectTableFileName, opts);
axialLength = subjectTable.Axial_Length_average;

% Load the anat measures
opts = detectImportOptions(p.Results.anatMeasuresFileName);
anatMeasuresTable = readtable(p.Results.anatMeasuresFileName, opts);

% The number of PCA coefficients we will retain
nCoeff = p.Results.nCoeff;
variableNames = {'AOSO_ID','PCA1_resid'};
for kk = 1:nCoeff
    variableNames = [variableNames {['PCA' num2str(kk)]}];
end

% This is the mmPerDeg at the ellipsoidal pole of the vitreous chamber
mmPerDeg = @(axialLength) (0.0165.*axialLength)-0.1070;

% Loop over layer sets
for ii = 1:length(p.Results.layerSetLabels)
    
    % Create a table to hold the PCA results
    pcaResultsTable = table('Size',[nSubs,nCoeff+2],'VariableTypes',repmat({'double'},1,nCoeff+2),'VariableNames',variableNames);
    
    % Loop over subjects and load the maps
    for ss = 1:nSubs
        fileName = fullfile(rawSubjectList(ss).folder,rawSubjectList(ss).name);
        load(fileName,'averageMaps');
        thisMap = averageMaps.(p.Results.layerSetLabels{ii});
        if ss==1
            imageSize = size(thisMap);
            avgMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
        end
        % Convert the map from thickness to volume per degree^2
        thisMap = thisMap .* mmPerDeg(axialLength(ss))^2;
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
    end
    
    % Create a 3D plot of the scores
    if p.Results.showPlots
        figure
        plot3(score(:,1),score(:,2),score(:,3),'*k')
        hold on
        plot3(score(subjectIdx(1),1),score(subjectIdx(1),2),score(subjectIdx(1),3),'*r')
        plot3(score(subjectIdx(2),1),score(subjectIdx(2),2),score(subjectIdx(2),3),'*b')
        axis equal
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
    end
    
    % Obtain the residuals of the first component after removing the effect
    % of axial length
    X = [ones(nSubs,1) axialLength-mean(axialLength)];
    B = X\score(:,1);
    B(1)=0;
    PCA1_resid = score(:,1) - X*B;
    pcaResultsTable.PCA1_resid = PCA1_resid;        

    % Now show images of the first, second, and third PCA components
    coeffMap = nan(imageSize(1),imageSize(2));
    
    if p.Results.showPlots
        figure
        for jj=1:nCoeff
            subplot(2,2,jj);
            coeffMap(observedIdx)=coeff(:,jj);
            mesh(coeffMap);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            axis square
            title('PCA coefficients (raw)')
        end
        figure
        for jj=1:nCoeff
            subplot(2,2,jj);
            coeffMap(observedIdx)=coeff(:,jj);
            
            
            scaleDown = 12;
            supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1)/scaleDown);
            downSampMap = coeffMap(rgcIplThicknessMap,1/scaleDown);
            [xo,yo,zo]=prepareSurfaceData(supportDeg,supportDeg,downSampMap);
            sf = fit([xo, yo],zo,'thinplateinterp');
            supportDeg = linspace(-p.Results.octRadialDegreesVisualExtent,p.Results.octRadialDegreesVisualExtent,imageSize(1));
            [xi,yi]=meshgrid(supportDeg,supportDeg);
            inRangeIdx = sqrt(xi.^2+yi.^2)<=octRadialDegreesVisualExtent;
            fitCoeffMap = sf(xi(inRangeIdx),yi(inRangeIdx));            
            mesh(fitCoeffMap);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            axis square
            title('PCA coefficients (raw)')
        end
    end
    
    % Show some PCA reconstructed OCT scans vs. the actual data
    if p.Results.showPlots
        figA = figure();
        figB = figure();
        zMax = 10;
        for jj = 1:length(subjectIdx)
            subjectID = split(rawSubjectList(jj).name,'_');
            subjectID = subjectID{1};
            
            figure(figA);
            subplot(length(subjectIdx),3,1+3*(jj-1))
            theMap = squeeze(avgMapBySubject(subjectIdx(jj),:,:));
            mesh(theMap);
            caxis([0 zMax]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([0 zMax]);
            title(['Subject ' subjectID] )
            axis square
            
            subplot(length(subjectIdx),3,2+3*(jj-1))
            reconMap = nan(imageSize(1),imageSize(2));
            reconMap(observedIdx)=Xfit(subjectIdx(jj),:);
            mesh(reconMap);
            caxis([0 zMax]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([0 zMax]);
            title('PCA reconstructed')
            axis square

            subplot(length(subjectIdx),3,3+3*(jj-1))
            mesh(reconMap-theMap);
            caxis([-zMax/2 zMax/2]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([-zMax/2 zMax/2]);
            title('Error')
            axis square

            figure(figB);
            plot(reconMap(imageSize(1)/2,:));
            hold on
            xlim([0 imageSize(2)]);
            ylim([0 zMax]);
            title('Horizontal meridian')
            axis square
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
