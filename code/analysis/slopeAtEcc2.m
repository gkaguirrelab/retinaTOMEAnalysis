function [YData, YData_errTop, YData_errBot]  = slopeAtEcc2(measurement,measurement2,normalizeFlag,interceptFlag)

 linearcoeff=zeros(2,size(measurement,1));
    linearcoeff_err=zeros(2,2,size(measurement,1));
    pvalue = zeros(1,size(measurement,1));
    for i = 1:size(measurement,1)
        [b,bint,r,rint,stats] = regress(measurement(i,:)', [measurement2(i,:)' ones(size(measurement2(i,:)'))]);
        linearcoeff(:,i)=b;
        linearcoeff_err(:,:,i)=bint;
        pvalue(i) = stats(3);
    end
    if(interceptFlag)
        s=2;
    else
        s=1;
    end
    err = squeeze(linearcoeff_err(s,2,:))-squeeze(linearcoeff(s,:))';
    ind = find(err'>10000|isnan(err')|err'==0);
    normalize = nanmean(measurement,2);
    %normalize = nanstd(measurement,0,2);

    if(~normalizeFlag)
        normalize = 1;
    end
    YData = linearcoeff(s,:)'./normalize;
    YData(ind) = nan;
    YData_errTop = squeeze(linearcoeff_err(s,1,:))./normalize;
    YData_errTop(ind) = nan;
    YData_errBot = squeeze(linearcoeff_err(s,2,:))./normalize;
    YData_errBot(ind) = nan;