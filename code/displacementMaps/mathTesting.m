clear all

%% 


%% Fit Functions
%input angle
angle = 0;
% fit the RGC density
[ecc_mm,RGCfit] = fitRGCdensity(angle,6);
% fit the RF density -- need to convert mm to deg
RFfit = fitRFdensity(convert_mm_to_deg(ecc_mm),angle);



%%
a =   1.809e+04
b =      -2.368
c =        8235
d =     -0.5072
syms ecc
RF_eq = @(ecc)(a*exp(b*ecc) + c*exp(d*ecc));

y_int = RF_eq(0);

RF_integ = int(RF_eq,ecc);

RF_integ = @(ecc)(- (1130625*exp(-(296*ecc)/125))/148 - (5146875*exp(-(317*ecc)/625))/317);

C = y_int - RF_integ(0); % 5.0201e+04

RF_integ = @(ecc)((- (1130625*exp(-(296*ecc)/125))/148 - (5146875*exp(-(317*ecc)/625))/317)+5.0201e+04);

% check if the integral and constant are correct
check = diff(RF_integ,ecc);

check = @(ecc)(18090*exp(-(296*ecc)/125) + 8235*exp(-(317*ecc)/625))



%%


[fitParams, fitRGCdensity] = fitFrechnetToRGCDensity(ecc_mm, RGCfit(ecc_mm), ones(size(ecc_mm))); 


RGC_int = @(ecc)((636320436822611.*ecc.^7)./64563604257983430656 - (39.*ecc.^6)./31250 - (1169.*ecc.^5)./500000 + ... 
          (3137.*ecc.^4)./40000 - (9649.*ecc.^3)./30000 - (3399.*ecc.^2)./10000 + (258.*ecc)./25 +10.32);
      
smps =log(0:0.5:20);      
v = exp(RGC_int(smps));
plot(0:0.5:20,v)



