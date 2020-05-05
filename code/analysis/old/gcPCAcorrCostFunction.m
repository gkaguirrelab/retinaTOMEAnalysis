function [slope, adjustedgcVolumePerDegSq] = gcPCAcorrCostFunction(k,gcVolumePerDegSq,gcVec,AreaPerDegSq,AxialLength,XPos_Degs)

%subtract by fraction of first PC for each subject
%  gcVec_zero = gcVec;
%  gcVec_zero(isnan(gcVec)) = 0;%need to zero nan's for PCA
%  
%  gcVec_new = zeros(size(gcVec));
%  
%  %fit a spline to remove non-zeros
%  for i = 1:size(gcVec,2)
%      x=1:size(gcVec,1);
%      currGC =  gcVec(:,i);
%      p = polyfit(x(~isnan(currGC))',currGC(~isnan(currGC)),15);
%      gcVec_new(:,i) = polyval(p,x');
%  end
%  
%  
% [coeff_thic,score_thic,~,~,explained_thic,mu_thic] = pca(gcVec_zero','Centered',false);
% %adjustedgcVolumePerDegSq = gcVolumePerDegSq' - (k*(coeff_thic(:,1)*score_thic(:,1)').*AreaPerDegSq)';
% adjustedgcVolumePerDegSq = gcVolumePerDegSq' - (k*(coeff_thic(:,1)*ones(50,1)').*AreaPerDegSq)';
% 
% PC1 = coeff_thic(:,1);

%subtract by fraction of each thickness profile
adjustedgcVolumePerDegSq = (gcVec.*(ones(size(gcVec)) - ones(size(gcVec,1),1)*(k./nanmean(gcVec,1))).*AreaPerDegSq)';
%Sanity check that above is equivalent to doing it line by line (Checked - it is)
% adjustedgcVolumePerDegSq2 = zeros(size(gcVec));
% for q=1:size(gcVec,2)
% adjustedgcVolumePerDegSq2(:,q) = (gcVec(:,q) - k*gcVec(:,q)/nanmean(gcVec(:,q))).*AreaPerDegSq(:,q);
% end
% adjustedgcVolumePerDegSq2=adjustedgcVolumePerDegSq2';
% max(max(abs(adjustedgcVolumePerDegSq-adjustedgcVolumePerDegSq2)))


%subtract by fraction of mean thickness profile
%adjustedgcVolumePerDegSq = gcVolumePerDegSq' - (k*(nanmean(gcVec,2)*ones(1,50)).*AreaPerDegSq)';
%Sanity Check: this is equivalent to above
%adjustedgcVolumePerDegSq = ((gcVec - k*(nanmean(gcVec,2)*ones(1,50))).*AreaPerDegSq)';


%subtract by a constant stromal thickness
%adjustedgcVolumePerDegSq = gcVolumePerDegSq' - k*AreaPerDegSq';

%PCA anlaysis of the volume
%adjustedgcVolumePerDegSq(isnan(gcVolumePerDegSq_zeroNaN))=0;
%[coeff,score,~,~,explained,mu] = pca(adjustedgcVolumePerDegSq,'Centered',true);
%[b,bint,r,rint,stats] = regress(score(:,1), [AxialLength ones(size(AxialLength))]);

%find median and it's slope with axial length
medianGC=nanmedian(adjustedgcVolumePerDegSq,2);
b=polyfit(AxialLength,medianGC,1);

%find median and it's slope with axial length
meanGC=nanmean(adjustedgcVolumePerDegSq,2);
b=polyfit(AxialLength,meanGC,1);
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