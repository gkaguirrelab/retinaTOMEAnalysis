function [slope, adjustedgcVolumePerDegSq] = gcPCAcorrCostFunction(k,gcVolumePerDegSq,gcVec,AreaPerDegSq,AxialLength,XPos_Degs)

%subtract by fraction of first PC for each subject
gcVec_zero = gcVec;
gcVec_zero(isnan(gcVec)) = 0;%need to zero nan's for PCA
[coeff_thic,score_thic,~,~,explained_thic,mu_thic] = pca(gcVec_zero','Centered',false);
adjustedgcVolumePerDegSq = gcVolumePerDegSq' - (k*(coeff_thic(:,1)*score_thic(:,1)').*AreaPerDegSq)';


%subtract by fraction of mean thickness profile
%adjustedgcVolumePerDegSq = gcVolumePerDegSq' - (k*(nanmean(gcVec,2)*ones(1,50)).*AreaPerDegSq)';

%subtract by a constant stromal thickness
%adjustedgcVolumePerDegSq = gcVolumePerDegSq' - k*AreaPerDegSq';

%PCA anlaysis of the volume
%adjustedgcVolumePerDegSq(isnan(gcVolumePerDegSq_zeroNaN))=0;
%[coeff,score,~,~,explained,mu] = pca(adjustedgcVolumePerDegSq,'Centered',true);
%[b,bint,r,rint,stats] = regress(score(:,1), [AxialLength ones(size(AxialLength))]);

%find median and it's slope with axial length
medianGC=nanmedian(adjustedgcVolumePerDegSq,2);
b=polyfit(AxialLength,medianGC,1);
slope=b(1);

%figure(30)
%hold on
% subplot(1,2,1)
% plot(k,slope);
% subplot(1,2,2)
%scatter(AxialLength,medianGC);
%plot(XPos_Degs,k*(nanmean(gcVec,2)))
%regressionPlot(AxialLength,medianGC, 'Axial length [mm]','Median GC Volume ', ...
%   [num2str(k) ' ' num2str(slope)],1,1)

%pause

end