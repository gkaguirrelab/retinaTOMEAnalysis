function density = coneDensityModel(x,angle,maxX,p)

a = p(1);
b = p(2);
c = p(3);
d = p(4);
ph1 = p(5);
f11 = p(6);
f12 = p(7);
f13 = p(8);
ph2 = p(9); % Unused params. The sin1 has the same angle as cos1
f21 = p(10);
f22 = p(11);
f23 = p(12);
ph3 = p(13);
f31 = p(14);
f32 = p(15);
f33 = p(16);
ph4 = p(17);
f41 = p(18);
f42 = p(19);
f43 = p(20);


g1 = f11.*gampdf(x,f12,f13)./max(gampdf(0:0.01:maxX,f12,f13));
g2 = f21.*gampdf(x,f22,f23)./max(gampdf(0:0.01:maxX,f22,f23));
g3 = f31.*gampdf(x,f32,f33)./max(gampdf(0:0.01:maxX,f32,f33));
g4 = f41.*gampdf(x,f42,f43)./max(gampdf(0:0.01:maxX,f42,f43));

m = g1.*cosd(angle+ph1) + ...
    g2.*sind(angle+ph1) + ...
    g3.*cosd(angle.*2+ph3) + ...
    g4.*cosd(angle.*4+ph4);

density = (m+1).*(a.*exp(b.*x)+c.*exp(d.*x));

end

