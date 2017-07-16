curcio_data = fullfile([getpref('octAnalysisForTOME','LocalDataPath') , '/Curcio/curcio_4meridian.mat']);
data=load(curcio_data);
data=data.data;


ecc_mm  = data(:,1);
ecc_deg = convert_mm_to_deg(ecc_mm);

rgcDensity_degSq_temporal = convert_mmSq_to_degSq(ecc_deg,data(:,2));
rgcDensity_degSq_superior = convert_mmSq_to_degSq(ecc_deg,data(:,4));
rgcDensity_degSq_nasal    = convert_mmSq_to_degSq(ecc_deg,data(:,6));
rgcDensity_degSq_inferior = convert_mmSq_to_degSq(ecc_deg,data(:,8));

f0 = 0.8928;
rm = 41.03;

midgetFractionByEccen = @(ecc_deg)(f0.*(1+(ecc_deg/rm)).^-1);

midgetFraction = midgetFractionByEccen(ecc_deg);
midget_rgcDensity_mmSq_temporal = rgcDensity_degSq_temporal.* midgetFraction;
midget_rgcDensity_mmSq_superior = rgcDensity_degSq_superior.* midgetFraction;
midget_rgcDensity_mmSq_nasal = rgcDensity_degSq_nasal.* midgetFraction;
midget_rgcDensity_mmSq_inferior = rgcDensity_degSq_inferior.* midgetFraction;

scaleData = 2*max([midget_rgcDensity_mmSq_temporal; midget_rgcDensity_mmSq_superior;...
    midget_rgcDensity_mmSq_nasal; midget_rgcDensity_mmSq_inferior]);

% norm_rgcDensity_temporal = midget_rgcDensity_mmSq_temporal./scaleData;
norm_rgcDensity_superior = midget_rgcDensity_mmSq_superior./scaleData;
% norm_rgcDensity_nasal = midget_rgcDensity_mmSq_nasal./scaleData;
% norm_rgcDensity_inferior = midget_rgcDensity_mmSq_inferior./scaleData;


[outParams, RGCdensityFit] = fitFrechetToRGCDensity(ecc_deg, norm_rgcDensity_superior, ones(size(ecc_deg)));



a    = outParams_RGC(1);
s    = outParams_RGC(2);
m    = outParams_RGC(3);



M = @(x)(2.*pi.*((x.*exp(-(((x+1)-m)/s).^-a) - (a.^-1).*(((x+1)-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,(((x+1)-m)/s).^a)))) - ...
    2.*pi.*((x.*exp(-((x-m)/s).^-a) - (a.^-1).*((x-m).*(((x-m)/s).^-a).^(-1/a).*gammainc(-1/a,((x-m)/s).^a)))));


figure
subplot(1,2,1)
plot(1:0.2:20,M(1:0.2:20).*scaleData)
title('RGC Cell Count')
xlabel('eccentricity (deg)')
ylabel('RGC Count (cells)')
subplot(1,2,2)
plot(1:0.2:20,cumsum(M(1:0.2:20).*scaleData))
title('Cumulative RGC Cell Count')
xlabel('eccentricity (deg)')
ylabel('Cumulative RGC Count (cells)')