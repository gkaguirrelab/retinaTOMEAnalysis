function [RFdensity,sampleBase_RF_deg] = densityRf(radDeg,smpPerDeg,verbose)
% Genrates a map of recetive feild dentity 
%
% From Eq. 7 in Drasdo et al. 2007 Vision Research.
% Formula is in Degrees 
% 
% The Drasdo equations are cast in terms of position in the visual field in
% degrees eccentricity and counts per square degree.
%
% The output of this routine will be in units of square degrees at an
% eccentricity given in degrees. The values, however, will be for that
% position in the retinal coordinate space, not the visual.
%
% The values returned here can be compared with Table 5 of Drasdo 2007,
% with the caveat that (e.g.) the superior meridian in the retina returned
% by this routine should be compared to the inferior meridian in the visual
% field presented in that table of Drasdo.


[~,theta,r] = createGrid(radDeg,smpPerDeg);

Rv0 = 0.011785;
Ro0 = 0.008333;

E20 = 20;
E2v = zeros(size(theta));
for i = 1:size(theta,1);
    for ii = 1:size(theta,2)
        if theta(i,ii) >= 0 & theta(i,ii) <= 90 
            E2v(i,ii) = interp1([0,90],[2.19,1.98],theta(i,ii));
         
        elseif theta(i,ii) > 90 & theta(i,ii) <= 180 
            E2v(i,ii) = interp1([90,180],[1.98,2.26],theta(i,ii));
          
        elseif theta(i,ii) > 180 & theta(i,ii) <= 270
            E2v(i,ii) = interp1([180,270],[2.26,1.5],theta(i,ii));
           
        elseif theta(i,ii) > 270 & theta(i,ii) <=360
            E2v(i,ii) = interp1([270,360],[1.5,2.19],theta(i,ii));
          
        end
    end
end
Rve = Rv0.*(1+r./E2v);
Roe =Ro0.*(1+r./E20);

k=1+(1.004 - 0.007209.*r + 0.001694.*r.^2 - 0.00003765.*r.^3).^-2;

RFdensity = ((1.12+0.0273.*r).*k)./(1.155.*(Rve.^2-Roe.^2));

%% Mask the output to not include vaules past the input radius 
mask = (r <= radDeg);
mask = double(mask);
mask(mask ==0) = nan;
RFdensity = RFdensity.*mask;

% Calculate a sample base
sampleBase_RF_deg = (0:1/smpPerDeg:radDeg); % the eccentricy of each sample of sup_RFdensity in degrees 

% This final "flip" of the image results in these relationships in the
% image:
%
%            midPoint = center pixel of Df
%            Temporal retina = left of midPoint of RFdensity
%            Nasal Retina    = right of midPoint of RFdensity
%            Supior Retina   = above midPoint of RFdensity
%            Inferior Retina = below midPoint of RFdensity
RFdensity=fliplr(RFdensity);

%% Validate the Output

if strcmp(verbose,'full')
    % 0-Nasal 90-Inferior 180-Temporal 270-Superior
    mdPt = round(size(RFdensity,1)/2);
    xSmpDegPos = 0:1/smpPerDeg:radDeg;
    xSmpDegNeg = radDeg:-1/smpPerDeg:0;
    figure;hold on;
    plot(xSmpDegNeg,RFdensity(mdPt,1:mdPt),'r')%Temporal
    plot(xSmpDegPos,RFdensity(mdPt:end,mdPt)','b')%Superior
    plot(xSmpDegPos,RFdensity(mdPt,mdPt:end),'g')%Nasal
    plot(xSmpDegNeg,RFdensity(1:mdPt,mdPt)','k')%Inferior
    legend('Temporal','Superior','Nasal','Inferior')
    xlabel('Eccentricity (deg)'); set(gca,'XScale','log');
    ylabel('Denstiy (deg^{-2}');  set(gca,'YScale','log');

end 
    
end