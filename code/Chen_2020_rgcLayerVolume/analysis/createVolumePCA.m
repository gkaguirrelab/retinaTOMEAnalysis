function [GCVolPCAScoreExpanded, GCVolPCAScoreExpandedSmoothed, GCVolPCACoeff, GCVolPCAVarExplained] = createVolumePCA(gcVolumePerDegSq,badIdx,XPos_Degs)

% Before we do the PCA, nan out the "bad" indices from the data. We won't
% attempt to reconstruct these points.
gcVolumePerDegSq(badIdx,:) = nan;

% We limit the PCA analysis to those x positions for which we have
% measurements for greater than 90% of the subjects.
nanX = sum(isnan(gcVolumePerDegSq'))>45;
gcVolumeCleaned = gcVolumePerDegSq(~nanX,:);

% Conduct a PCA using the alternating least squares (ALS) algorithm to
% account for the few missing values
[GCVolPCACoeff,score,~,~,GCVolPCAVarExplained] = pca(gcVolumeCleaned,'Centered',false,'algorithm','als');

% Expand the score vectors back to the full x range
GCVolPCAScoreExpanded = nan(size(gcVolumePerDegSq));
GCVolPCAScoreExpanded(~nanX,:) = score;

% Find the three segment domains
tt = find(diff(nanX));
domains = {[1 tt(1)],[tt(2)+1 tt(3)],[tt(4)+1 length(nanX)]};

% Perform piece-wise spline smoothing of the scores to remove the noisy
% effects of data imputation
smoothVal = 0.1; % 0-1, lower is smoother.
GCVolPCAScoreExpandedSmoothed = GCVolPCAScoreExpanded;
for cc = 1:size(GCVolPCAScoreExpanded,2)
    for dd = 1:length(domains)
        rd = domains{dd};
        x = XPos_Degs(rd(1):rd(2));
        y = GCVolPCAScoreExpanded(rd(1):rd(2),cc);
        pp = csaps(x',y,smoothVal);
        GCVolPCAScoreExpandedSmoothed(rd(1):rd(2),cc) = ppval(pp,x');
    end
end


