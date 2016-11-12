%this processes the curcio RGCd 4 meridian data

load('curcio_4meridian.mat')

ecc_mm = data(:,1);
temp = data(:,2);
sup = data(:,4);
nasal = data(:,6);
inferior = data(:,8);

% Conversion of eccentricities in millimeters to degrees; Watson 2014

ecc_deg = 3.556*ecc_mm+0.05992*ecc_mm.^2-0.007358*ecc_mm.^3+0.0003027*ecc_mm.^4;



alpha= 0.0752+5.846e-5*ecc_mm-1.064e-5*ecc_mm.^2+4.116e-8*ecc_mm.^3;

ecc_thresh = 23;


temp_deg2 = temp.*alpha;
sup_deg2 = sup.*alpha;
nasal_deg2 = nasal.*alpha;
inferior_deg2 = inferior.*alpha;


%%

[curve_nasal] = fit(ecc_deg,nasal_deg2 ,'smoothingspline','Exclude', find(isnan(nasal_deg2)),'SmoothingParam', 1);

[curve_sup] = fit(ecc_deg,sup_deg2,'smoothingspline', 'Exclude',find(isnan(sup_deg2)),'SmoothingParam', 1);

[curve_temp] = fit(ecc_deg,temp_deg2,'smoothingspline', 'Exclude',find(isnan(temp_deg2)),'SmoothingParam', 1);

[curve_inferior] = fit(ecc_deg,inferior_deg2,'smoothingspline', 'Exclude',find(isnan(inferior_deg2)),'SmoothingParam', 1);
%%Check the Spline
figure
hold on 
plot(ecc_deg,curve_temp(ecc_deg),'r')%temp
scatter(ecc_deg,temp_deg2,'r')
plot(ecc_deg,curve_sup(ecc_deg),'b'); %sup
scatter(ecc_deg,sup_deg2,'b')
plot(ecc_deg,curve_nasal(ecc_deg),'g'); %nasal
scatter(ecc_deg,nasal_deg2,'g')
plot(ecc_deg,curve_inferior(ecc_deg),'k'); %inf
scatter(ecc_deg,inferior_deg2,'k')
set(gca,'xscale','log')
set(gca,'yscale','log')



%%

eccSmp = 20;
bkgrd = zeros(2*eccSmp+1,2*eccSmp+1);
temp_hm = curve_temp(0:eccSmp);
sup_vm = curve_sup(eccSmp:-1:0)';
nasal_hm = curve_nasal(eccSmp:-1:0);
inf_vm = curve_inferior(0:eccSmp)';

midPt = round((2*eccSmp)/2)+1;
bkgrd(midPt,midPt:end) = temp_hm;
bkgrd(midPt,1:midPt) = nasal_hm;
bkgrd(midPt:end,midPt) = inf_vm;
bkgrd(1:midPt,midPt) = sup_vm;



%% half size shift matrix to zero center the image


