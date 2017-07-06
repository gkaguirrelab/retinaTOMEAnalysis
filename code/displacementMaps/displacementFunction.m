function Xc = displacementFunction(Xr,params) 
a = params.a;
b = params.b;
c = params.c;
d = params.d;

s = params.s;
m = params.m;
alpha = params.alpha;
g = params.g;

Kr = params.Kr;
Kc = params.Kc;


Func = log((a.*exp(-b.*Xr)./(g.*b))+(c.*exp(-d.*Xr)./(g.*d))+((Kr+Kc)/g)).^-1;

Xc = -s.*nthroot(Func,alpha)-m;

end