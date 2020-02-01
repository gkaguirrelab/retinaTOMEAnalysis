function slope = gcPCAcorrCostFunction(k,gcVolumePerDegSq,gcVec,AreaPerDegSq,AxialLength)


gcVolumePerDegSq_zeroNaN = gcVolumePerDegSq' - (k*(nanmean(gcVec,2)*ones(1,50)).*AreaPerDegSq)';

%gcVolumePerDegSq_zeroNaN = gcVolumePerDegSq' - k*AreaPerDegSq';
%gcVolumePerDegSq_zeroNaN = gcVolumePerDegSq' - (k*nanmean(gcVolumePerDegSq,2)*ones(1,50))';
%gcVolumePerDegSq_zeroNaN(isnan(gcVolumePerDegSq_zeroNaN))=0;

%[coeff,score,~,~,explained,mu] = pca(gcVolumePerDegSq_zeroNaN,'Centered',true);
%[b,bint,r,rint,stats] = regress(score(:,1), [AxialLength ones(size(AxialLength))]);

medianGC=nanmedian(gcVolumePerDegSq_zeroNaN,2);
%[b,bint,r,rint,stats] = regress(medianGC, [AxialLength ones(size(AxialLength))]);
b=polyfit(AxialLength,medianGC,1);

% meanGC=nanmean(gcVolumePerDegSq_zeroNaN,2);
% %[b,bint,r,rint,stats] = regress(medianGC, [AxialLength ones(size(AxialLength))]);
% b=polyfit(AxialLength,meanGC,1);


slope=b(1);

figure(2)
% subplot(1,2,1)
% plot(k,slope);
% subplot(1,2,2)
%scatter(AxialLength,medianGC);

%regressionPlot(AxialLength,medianGC, 'Axial length [mm]','Median GC Volume ', ...
%    [num2str(k) ' ' num2str(slope)],1)

%pause

end