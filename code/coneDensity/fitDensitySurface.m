

% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')



% Load the confocal and split data
splitFiles=dir('densityAnalysis/*_split.mat');
splitNames = strrep(extractfield(splitFiles,'name'),'_split.mat','');

confocalFiles=dir('densityAnalysis/*_confocal.mat');
confocalNames = strrep(extractfield(confocalFiles,'name'),'_confocal.mat','');

supportLength = 1799;
maxSupportDeg = 15;
supportDegDelta = 0.0078;
conStart = 0.5;
conStop = 1.6;
splitStart = 2.1;
nFourier = 3;


% Define the eccentricity support, and the ranges (in degrees) that will be
% used for the confocal and split detecton data sets
supportDeg = 0:supportDegDelta:supportDegDelta*(supportLength-1);
idxA = find(supportDeg>=conStart,1);
idxB = find(supportDeg>=conStop,1);
idxC = find(supportDeg>=splitStart,1);


%% Loop through subjects and create the composite polar density image

subNames = unique([splitNames confocalNames]);
dataMat = nan(supportLength,supportLength,length(subNames));

for ss = 1:length(subNames)
    
    % A matrix to hold the data for this subject
    y = nan(supportLength,supportLength);
    
    % Add the confocal data
    conFile = fullfile('densityAnalysis',[subNames{ss} '_confocal.mat']);
    if isfile(conFile)
        load(conFile,'data');
        y(:,idxA:idxB) = data.polarDensity(:,idxA:idxB);
    end
    
    % Add the split data
    splitFile = fullfile('densityAnalysis',[subNames{ss} '_split.mat']);
    if isfile(splitFile)
        load(splitFile,'data');
        y(:,idxC:end) = data.polarDensity(:,idxC:end);
    end
    
    % Filter out any negative values
    y(y<0)=nan;
    dataMat(:,:,ss)=y;
    
end

polarRatio = (size(dataMat,1)+1)/360;

%% Initial fit to mean density profile

% Density, eccentricity, polar angle
Y = nanmean(nanmean(dataMat,3));
w = nanmean(sum(~isnan(dataMat),3));
X = supportDeg;
P = zeros(size(X));

% objective
validIdx = ~isnan(Y);
myObj = @(p) norm( w(validIdx).* (Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );

% p0 and bounds
asymptote = 0; % In case we want to force a minimum asymptote value
pBlock0 = [1300, -0.25, 8500, -1.25, asymptote];
pBlockLB = [0,-5,0,-5,asymptote];
pBlockUB = [5e4,0,5e4,0,asymptote];
mBlock0 = [0 0 1.01 1.01];

p0 = [pBlock0, repmat(mBlock0,1,nFourier)];
lb = [pBlockLB, repmat(mBlock0,1,nFourier)];
ub = [pBlockUB, repmat(mBlock0,1,nFourier)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub);


%% Fit polar angle variation

% Density, eccentricity, polar angle
Y = nanmean(dataMat,3);
w = sum(~isnan(dataMat),3);
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);

% objective
validIdx = ~isnan(Y);
myObj = @(p) norm( w(validIdx).* (Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );

% p0 and bounds
mBlockLB = [-45 -1 2 0.01];
mBlockUB = [45 1 10 4];

p0 = [p(1:5), repmat(mBlock0,1,nFourier)];
lb = [pBlockLB, repmat(mBlockLB,1,nFourier)];
ub = [pBlockUB, repmat(mBlockUB,1,nFourier)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub);


%% Fit an individual subject

% % Data
% Y = squeeze(dataMat(:,:,3));
%
% % objective
% validIdx = ~isnan(Y);
% myObj = @(p) norm( Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p) );
%
% % p0 and bounds
% mBlockFit = p(6:end);
% p0 = [p(1:5), mBlockFit];
% lb = [pBlockLB, mBlockFit];
% ub = [pBlockUB, mBlockFit];
%
% % search
% pSub1 = fmincon(myObj,p0,[],[],[],[],lb,ub);

foo=1;



%% plot

% generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:)=myModel(X,P,maxSupportDeg,p);

meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];

figure
surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none');
hold on
plot3(X(:),P(:),Y(:),'.k')
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');


figure
for ii = [0.375 0.75 1.5 3 6 10]
    idx = find(supportDeg>ii,1);
    semilogy(Y(:,idx),'.');
    hold on
    semilogy(Yfit(:,idx),'-r');
    text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
end
xticks(meridianAngles*polarRatio);
xticklabels(meridianLabels);
ylim([10^2,10^4]);
ylabel('log_1_0 density [cones/deg^2]')


figure
for mm=1:4
    subplot(2,2,mm)
    plot(supportDeg,Y(round((meridianAngles(mm))*polarRatio+1),:),'.k');
    hold on
    plot(supportDeg,Yfit(round((meridianAngles(mm))*polarRatio+1),:),'-r');
    xlabel('Eccentricity [deg]');
    ylabel('Density [cones/deg^2]');
    title(meridianLabels{mm});
end

figure
meridianSpec = {'-r','--b','-m','--g'};
for mm=1:4
    plot(supportDeg,Yfit(meridianAngles(mm)*polarRatio+1,:),meridianSpec{mm},'LineWidth',2);
    hold on
end
xlabel('Eccentricity [deg]');
ylabel('Density [cones/deg^2]');
legend(meridianLabels(1:4));


figure
for gg = 1:nFourier
    ph = p((gg-1)*4+6);
    f1 = p((gg-1)*4+7);
    f2 = p((gg-1)*4+8);
    f3 = p((gg-1)*4+9);
    g = f1.*gampdf(supportDeg,f2,f3)./max(gampdf(0:0.01:maxSupportDeg,f2,f3));
    plot(supportDeg,g)
    hold on
end
lengendLabels = {'cos1','cos2','cos4'};
legend(lengendLabels(1:nFourier));
xlabel('Eccentricity [deg]');
ylabel('Modulation [a.u.]');


figure
imagesc(w)
axis square
yticks(meridianAngles*polarRatio);
yticklabels(meridianLabels);
deltaDeg = supportDeg(2)-supportDeg(1);
xticks(round(0:1/deltaDeg:floor(maxSupportDeg)/deltaDeg));
xticklabels(0:1:floor(maxSupportDeg));
xlabel('Eccentricity [deg]');
colorbar
title('Weight map');


foo=1;



function density = myModel(x,angle,maxX,p)

a = p(1);
b = p(2);
c = p(3);
d = p(4);
e = p(5);
ph1 = p(6);
f11 = p(7);
f12 = p(8);
f13 = p(9);
ph2 = p(10);
f21 = p(11);
f22 = p(12);
f23 = p(13);
ph3 = p(14);
f31 = p(15);
f32 = p(16);
f33 = p(17);

g1 = f11.*gampdf(x,f12,f13)./max(gampdf(0:0.01:maxX,f12,f13));
g2 = f21.*gampdf(x,f22,f23)./max(gampdf(0:0.01:maxX,f22,f23));
g3 = f31.*gampdf(x,f32,f33)./max(gampdf(0:0.01:maxX,f32,f33));

m = g1.*cosd(angle+ph1) + ...
    g2.*cosd(angle.*2+ph2) + ...
    g3.*cosd(angle.*4+ph3);

density = (m+1).*(a.*exp(b.*x)+c.*exp(d.*x)+e);

end

