function [RFdensity,sampleBase_RF_deg] = densityRF_watson_eq(radDeg,smpPerDeg,interp,verbose)

sampleBase_RF_deg = 0:1/smpPerDeg:radDeg; 

nasal          = 2*(14804.6) .* (0.9729*((1+sampleBase_RF_deg./1.084)).^-2)+(1-0.9729).*exp(-1.*sampleBase_RF_deg./7.633);
nasal          = nasal .* (0.8928*((1+sampleBase_RF_deg./41.03).^-1));
curve_nasal    = fit(sampleBase_RF_deg',nasal','smoothingspline','Exclude', find(isnan(nasal)),'SmoothingParam', 1);

superior       = 2*(14804.6) * ( 0.9935*(1+sampleBase_RF_deg/(1.035)).^-2+(1-0.9935)*exp(-1*sampleBase_RF_deg/16.35));
superior          = superior .* (0.8928*(1+sampleBase_RF_deg./41.03).^-1);
curve_superior = fit(sampleBase_RF_deg',superior','smoothingspline','Exclude', find(isnan(superior)),'SmoothingParam', 1);

temporal       = 2*(14804.6) * ( 0.9851*(1+sampleBase_RF_deg/(1.058)).^-2+(1-0.9851)*exp(-1*sampleBase_RF_deg/22.14));
temporal          = temporal .* (0.8928*(1+sampleBase_RF_deg./41.03).^-1);
curve_temporal = fit(sampleBase_RF_deg',temporal','smoothingspline','Exclude', find(isnan(temporal)),'SmoothingParam', 1);

inferior       = 2*(14804.6) * ( 0.996*(1+sampleBase_RF_deg/(0.9932)).^-2+(1-0.996)*exp(-1*sampleBase_RF_deg/12.13));
inferior          = inferior .* (0.8928*(1+sampleBase_RF_deg./41.03).^-1);
curve_inferior = fit(sampleBase_RF_deg',inferior','smoothingspline','Exclude', find(isnan(inferior)),'SmoothingParam', 1);

[~,meridian(:,:,2),meridian(:,:,3)] = createGrid(radDeg,smpPerDeg);


%%interpolate 
for i = 1:size(meridian,1);
    for ii = 1:size(meridian,2)
        
        if meridian(i,ii,2) >= 0 & meridian(i,ii,2) <= 90
            VMd = curve_superior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_nasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([0,90],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) > 90 & meridian(i,ii,2) <= 180
            VMd = curve_superior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_temporal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([90,180],[VMd,HMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 180 & meridian(i,ii,2) <= 270
            VMd = curve_inferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_temporal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([180,270],[HMd,VMd],theta,interp);
            
        elseif meridian(i,ii,2) >= 270 & meridian(i,ii,2) <=360
            VMd = curve_inferior(meridian(i,ii,3));
            VMd(VMd<0) = 0;
            HMd = curve_nasal(meridian(i,ii,3));
            HMd(HMd<0) = 0;
            theta = meridian(i,ii,2);
            meridian(i,ii,1) = interp1([270,360],[VMd,HMd],theta,interp);
            
        end
    end
end 

mask=(meridian(:,:,3)<=radDeg);
double(mask(mask == 0)) = nan;
RFdensity =meridian(:,:,1).*mask;

%% Validate the Output

if strcmp(verbose,'full')
    % 0-Nasal 90-Surperior 180-Temporal 270-Inferior
    rgc = meridian(:,:,1);
    mdPt = round(size(rgc,1)/2);
    xSmpDegPos = 0:1/smpPerDeg:radDeg;
    xSmpDegNeg = radDeg:-1/smpPerDeg:0;
    figure;hold on;
    plot(xSmpDegNeg,rgc(mdPt,1:mdPt),'r')%Temporal
    plot(xSmpDegPos,rgc(mdPt:end,mdPt)','b')%Inferior
    plot(xSmpDegPos,rgc(mdPt,mdPt:end),'g')%Nasal
    plot(xSmpDegNeg,rgc(1:mdPt,mdPt)','k')%Superior 
    legend('Temporal','Inferior','Nasal','Superior')
    xlabel('Eccentricity (deg)'); set(gca,'XScale','log');
    ylabel('Denstiy (deg^{-2}');  set(gca,'YScale','log');

end 


end

