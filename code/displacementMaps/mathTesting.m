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
clear all

p1 =   6.899e-05;
p2 =   -0.007488;
p3 =    -0.01169;
p4 =      0.3137;
p5 =     -0.9649;
p6 =     -0.6798;
p7 =       10.32;

syms ecc
RGC_eq = @(ecc)(p1*ecc.^6 + p2*ecc.^5 + p3*ecc.^4 + p4*ecc.^3 + p5*ecc.^2 + p6*ecc + p7);

RGC_int = int(RGC_eq,ecc);


RGC_int = @(ecc)((636320436822611.*ecc.^7)./64563604257983430656 - (39.*ecc.^6)./31250 - (1169.*ecc.^5)./500000 + ... 
          (3137.*ecc.^4)./40000 - (9649.*ecc.^3)./30000 - (3399.*ecc.^2)./10000 + (258.*ecc)./25 +10.32);
      
smps =log(0:0.5:20);      
v = exp(RGC_int(smps));
plot(0:0.5:20,v)


% 
% 
% 
% %%
% a =   6.899e-05
% b =   -0.007488
% c =    -0.01169
% d =      0.3137
% e =     -0.9649
% f =     -0.6798
% g =       10.32
% 
% syms x
% RGC_eq = @(ecc)(p1*ecc.^6 + p2*ecc.^5 + p3*ecc.^4 + p4*ecc.^3 + p5*ecc.^2 + p6*ecc + p7);
% RGC_eq = @(x)exp((exp(h) .* (x.^6).^a .* (x.^5).^b .* (x.^4).^c  .* (x.^3).^d .* (x.^2).^f .* x.^g));
% 
% 
% y_value = RGC_eq(1)
% 
% 
% RGC_integ= @(x)(exp(g).*x.^(6.*a+5.*b+4.*c+3.*d+2.*e+f+1)./(6.*a+5.*b+4.*c+3.*d+2.*e+f+1));
% 
%             
% RGC_integ= @(x)(exp(-40316.1./x.^0.752386))             
%             
% C = y_value - RF_integ(1); % 5.0201e+04
% 
% RF_integ = @(ecc)((- (1130625*exp(-(296*ecc)/125))/148 - (5146875*exp(-(317*ecc)/625))/317)+5.0201e+04);
% 
% % check if the integral and constant are correct
% check = diff(RF_integ,ecc);
% 
% check = @(ecc)(18090*exp(-(296*ecc)/125) + 8235*exp(-(317*ecc)/625))
% 


