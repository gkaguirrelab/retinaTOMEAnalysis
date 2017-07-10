function Xc = displacementFunction(X_RF,params) 
a = params.a;
b = params.b;
c = params.c;
d = params.d;

s = params.s;
m = params.m;
alpha = params.alpha;
g = params.g;

K_RF = params.Kr;
K_RGC = params.Kc;

y = a.*exp(-b.*X_RF)+c.*exp(-d.*X_RF);

density_Term1 = (a.*exp(-1*b.*X_RF))./b;
density_Term2 = (c.*exp(-1*d.*X_RF))./d;
Density_RF = -1.*density_Term1 - density_Term1 + K_RF;

y = g*((alpha./s)*((X_RF-m)./s).^(-1-alpha).* exp( -((X_RF-m)./s).^(-alpha)));

A = log((Density_RF - K_RGC)./g).^-1;

-S.*nthroot(A,alpha)-m

Xrgc = 



end