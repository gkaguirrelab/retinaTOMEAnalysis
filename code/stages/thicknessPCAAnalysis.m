function thicknessPCAAnalysis(dataDir, varargin)
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

%% Parse and check the parameters
p.parse(dataDir, varargin{:});


% Obtain a list of subjects
rawSubjectList = dir(fullfile(dataDir,'*/*.mat'));
nSubs = length(rawSubjectList);

% Loop over layer sets
for ii = 1:length(p.Results.layerSetLabels)
    
    % Loop over subjects and load the maps
    for ss = 1:length(rawSubjectList)
        fileName = fullfile(rawSubjectList(ss).folder,rawSubjectList(ss).name);
        load(fileName,'averageMaps');
        thisMap = averageMaps.(p.Results.layerSetLabels{ii});
        if ss==1
            imageSize = size(thisMap);
            avgMapBySubject = nan(nSubs,imageSize(1),imageSize(2));
        end
        avgMapBySubject(ss,:,:)=thisMap;
    end
    
    % Get the mean thickness
    meanThicknessBySubject = squeeze(nanmean(nanmean(avgMapBySubject,2),3));
    
    % Determine which points are not nans
    observedIdx = ~isnan(squeeze(sum(avgMapBySubject,1)));
    
    % Conduct a PCA analysis across subjects, averaged over eyes. First,
    % Convert the non-nan portion of the image to a vector
    X = [];
    for ss=1:nSubs
        tmp = squeeze(avgMapBySubject(ss,:,:));
        X(ss,:) = tmp(observedIdx);
    end
    
    % Calc the PCA
    [coeff,score,~,~,explained,mu] = pca(X,'Centered',false);
    
    % Decide how many coefficients to keep
    explained = [nan; explained(2:end)./sum(explained(2:end))];
    nCoeff = 9;
    
    % Obtain the fit to the data
    Xfit = score(:,1:nCoeff)*coeff(:,1:nCoeff)'+(repmat(mu,nSubs,1));
    
    % Find the two subjects with the values close to the median for the
    % first coefficient, but maximally different from each other on the
    % second coefficient
    medianThickness = median(meanThicknessBySubject);
    [~,idxThickOrder] = sort(abs(meanThicknessBySubject-medianThickness),'ascend');
    rankThickOrder(idxThickOrder)=1:length(idxThickOrder);
    [~,idxLowCoeff] = sort(score(:,2),'ascend');
    rankLowCoeff(idxLowCoeff)=1:length(idxLowCoeff);
    [~,idxHighCoeff] = sort(score(:,2),'descend');
    rankHighCoeff(idxHighCoeff)=1:length(idxHighCoeff);
    
    rankThickLow = mean([rankThickOrder; rankLowCoeff]);
    rankThickHigh = mean([rankThickOrder; rankHighCoeff]);
    [~,subjectIdx(1)] = min(rankThickLow);
    [~,subjectIdx(2)] = min(rankThickHigh);
    
    
    % Show a plot of the explanatory power of the components
    if p.Results.showPlots
        figure
        plot(explained)
    end
    
    % Create a 3D plot of the scores for the first 3 components
    if p.Results.showPlots
        figure
        plot3(score(:,1),score(:,2),score(:,3),'*k')
        hold on
        plot3(score(subjectIdx(1),1),score(subjectIdx(1),2),score(subjectIdx(1),3),'*b')
        plot3(score(subjectIdx(2),1),score(subjectIdx(2),2),score(subjectIdx(2),3),'*r')
        axis equal
    end
    
    % Now show images of the first, second, and third PCA components
    coeffMap = nan(imageSize(1),imageSize(2));
    
    if p.Results.showPlots
        figure
        for jj=1:nCoeff
            subplot(3,3,jj);
            coeffMap(observedIdx)=coeff(:,jj);
            mesh(coeffMap);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            axis square
        end
    end
    
    % Show some PCA reconstructed OCT scans vs. the actual data
    if p.Results.showPlots
        for jj = 1:length(subjectIdx)
            subjectID = split(rawSubjectList(1).name,'_');
            subjectID = subjectID{1};
            
            subplot(length(subjectIdx),3,1+3*(jj-1))
            mesh(squeeze(avgMapBySubject(subjectIdx(jj),:,:)));
            caxis([0 125]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([0 150]);
            title(['Subject ' subjectID] )
            axis square
            
            subplot(length(subjectIdx),3,2+3*(jj-1))
            reconMap = nan(imageSize(1),imageSize(2));
            reconMap(observedIdx)=Xfit(subjectIdx(jj),:);
            mesh(reconMap);
            caxis([0 125]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([0 150]);
            title('PCA reconstructed')
            axis square
            
            subplot(length(subjectIdx),3,3+3*(jj-1))
            mesh(reconMap-squeeze(avgMapBySubject(subjectIdx(jj),:,:)));
            caxis([-62.5 62.5]);
            xlim([0 imageSize(1)]);
            ylim([0 imageSize(1)]);
            zlim([-75 75]);
            title('reconstruction error')
            axis square
        end
    end
    
end
