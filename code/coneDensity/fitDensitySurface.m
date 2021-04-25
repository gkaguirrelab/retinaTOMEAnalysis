function [p, Yfit, fVal] = fitDensitySurface(Y,w,preFitAvgEccen,simplePolarModel,p0,supportDeg,maxSupportDeg)
% Fit a multi-parameter surface to cone density data
%
% Syntax:
%   output = myFunc(input)
%
% Description:
%   Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean euismod
%   nulla a tempor scelerisque. Maecenas et lobortis est. Donec et turpis
%   sem. Sed fringilla in metus ut malesuada. Pellentesque nec eros
%   efficitur, pellentesque nisl vel, dapibus felis. Morbi eu gravida enim.
%   Sed sodales ipsum eget finibus dapibus. Fusce sagittis felis id orci
%   egestas, non convallis neque porttitor. Proin ut mi augue. Cras posuere
%   diam at purus dignissim, vel vestibulum tellus ultrices
%
% Inputs:
%   none
%   foo                   - Scalar. Foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo
%
% Optional key/value pairs:
%   none
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
%}


arguments
    Y (:,:) {mustBeNumeric}
    w (:,:) {mustBeNumeric} = ones(size(Y))
    preFitAvgEccen (1,1) = true;
    simplePolarModel (1,1) = true;
    p0 (1,20) {mustBeNumeric} = [1.4219e+03, -0.0813, 9.4727e+03, -1.2994, 18.9946, -0.0663, 10.1278, 0.1424, 0, -0.0627, 6.8191, 0.8296, -10.8437, 0.0925, 3.9524, 0.4707, -3.6584, 0.1173, 2.9143, 1.999]
    supportDeg (1,:) {mustBeNumeric} = 0:0.0078:0.0078*(size(Y,1)-1)
    maxSupportDeg (1,1) {mustBeNumeric} = 15
end

%% pBlock and mBlock settings
% These are the parameters of the model. The pBlock defines the form of the
% sum of exponential fit. The mBlock defines the parameters of a gamma
% function modulation of a sinusoidal variation in density across polar
% angle.

pBlockLB = [0,-5,0,-5];
pBlockUB = [5e4,0,5e4,0];

mBlockLB = [-35 -1 2 0.01];
mBlockUB = [35 1 20 2];

% The number of Fourier components that models variation across polar angle
nFourier = 4;

% The length of the support, pulled out here for code clarity
supportLength = length(supportDeg);


%% Initial fit to mean density profile

if preFitAvgEccen

    % Density, eccentricity, polar angle
    Ym = nanmean(Y,1);
    wm = nanmean(w,1);
    X = supportDeg;
    P = zeros(size(X))+90;
    
    % Bounds
    lb = [pBlockLB, p0(5:end)];
    ub = [pBlockUB, p0(5:end)];

    % Options
    options = optimoptions('fmincon','Display','off','Algorithm','interior-point');

    % objective
    validIdx = ~isnan(Ym);
    myObj = @(p) norm( wm(validIdx).* (Ym(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );
    
    % search
    pm = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % Update p0
    p0(1:4) = pm(1:4);

end


%% Fit polar angle variation

% Density, eccentricity, polar angle
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);
validIdx = ~isnan(Y);

if simplePolarModel
    
    % p0 and bounds
    lb = [pBlockLB, -30 0];
    ub = [pBlockUB, 30 2];
    
    % assemble p function
    pFull = @(p) [p(1:4) p0(5)+p(5) p0(6)*p(6) p0(7:8) p0(9)+p(5) p0(10)*p(6) p0(11:12) p0(13)+p(5) p0(14)*p(6) p0(15:16) p0(17)+p(5) p0(18)*p(6) p0(19:20)];
    
    % truncate p0
    p0 = [p0(1:4) 0 1];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,pFull(p))) );
    
else
    
    % p0 and bounds
    lb = [pBlockLB, repmat(mBlockLB,1,nFourier)];
    ub = [pBlockUB, repmat(mBlockUB,1,nFourier)];
    
    % objective
    myObj = @(p) norm( w(validIdx).* (Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );
    
end

% search
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','UseParallel',true);
[p, fVal] = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

if simplePolarModel
    p = pFull(p);
end

% generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:)=myModel(X,P,maxSupportDeg,p);


end


%% Local functions

function density = myModel(x,angle,maxX,p)

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

