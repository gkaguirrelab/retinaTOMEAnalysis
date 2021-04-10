

% The overal result directory
cd('/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/Connectome_AOmontages_images')

% Load the confocal and split data
dataFileName = 'confocalDensityProfileData.mat';
load(dataFileName)
confocalData = data;

dataFileName = 'splitDensityProfileData.mat';
load(dataFileName)
splitData = data;

clear data

conStart = 0.5;
conStop = 1.6;
splitStart = 2.1; 


%% Aggregate the profiles
% Find the longest support deg
supportLength = size(confocalData{1}.polarDensity,1);
supportDeg = 0:confocalData{1}.meta.supportDegDelta:confocalData{1}.meta.supportDegDelta*(supportLength-1);

idxA = find(supportDeg>=conStart,1);
idxB = find(supportDeg>=conStop,1);
idxC = find(supportDeg>=splitStart,1);

%% Loop through subjects and create the composite polar density image
subNames = unique(cellfun(@(x) x.meta.subName,[confocalData splitData],'UniformOutput',false));
dataMat = nan(supportLength,supportLength,length(subNames));
for ss = 1:length(subNames)
    y = nan(supportLength,supportLength);
    conIdx = find(cellfun(@(x) strcmp(x.meta.subName,subNames{ss}),confocalData));
    if ~isempty(conIdx)
        y(:,idxA:idxB) = confocalData{conIdx}.polarDensity(:,idxA:idxB);
    end
    splitIdx = find(cellfun(@(x) strcmp(x.meta.subName,subNames{ss}),splitData));
    if ~isempty(splitIdx)
        y(:,idxC:end) = splitData{splitIdx}.polarDensity(:,idxC:end);
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
maxSupportDeg = supportDeg(find(~isnan(Y), 1, 'last' ));
validIdx = ~isnan(Y);
myObj = @(p) norm( w(validIdx).* (Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p)) );

% p0 and bounds
asymptoteIdx = find(supportDeg>10,1);
asymptote = nanmean(Y(asymptoteIdx:end));


mBlock0 = [0 0 1.01 1.01];

p0 = [1300, -0.25, 8500, -1.25, asymptote, repmat(mBlock0,1,4)];
lb = [0,-5,0,-5,asymptote, repmat(mBlock0,1,4)];
ub = [5e4,0,5e4,0,asymptote, repmat(mBlock0,1,4)];

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
mBlockLB = [-180 -1 1.01 1.01];
mBlockUB = [180 1 20 20];

p0 = [p(1:5), repmat(mBlock0,1,4)];
lb = [p(1:5), repmat(mBlockLB,1,4)];
ub = [p(1:5), repmat(mBlockUB,1,4)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub);

% generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:)=myModel(X,P,maxSupportDeg,p);

%% plot

meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];

figure
surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none')
hold on
plot3(X(:),P(:),Y(:),'.k')
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');

figure
for ii = [0.375 0.75 1.5 3 6 12]
    idx = find(supportDeg>ii,1);
    semilogy(Y(:,idx),'.');
    hold on
    semilogy(Yfit(:,idx),'-r');
    text(300*polarRatio,Yfit(45*polarRatio,idx),sprintf('%2.1f°',supportDeg(idx)));
end
xticks(meridianAngles*polarRatio);
xticklabels(meridianLabels);
ylabel('log_1_0 density [cones/deg^2]')

figure
for mm=1:4
    subplot(2,2,mm)
    plot(supportDeg,Y(meridianAngles(mm)*polarRatio+1,:),'.k');
    hold on
    plot(supportDeg,Yfit(meridianAngles(mm)*polarRatio+1,:),'-r');  
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
for gg = 1:4
    ph = p((gg-1)*4+6);
    f1 = p((gg-1)*4+7);
    f2 = p((gg-1)*4+8);
    f3 = p((gg-1)*4+9);
    g = f1.*gampdf(supportDeg,f2,f3)./max(gampdf(0:0.01:maxSupportDeg,f2,f3));
    plot(supportDeg,g)
    hold on
end
legend({'sin1','sin2','sin4','sin8'});
xlabel('Eccentricity [deg]');
ylabel('Modulation [a.u.]');


figure
imagesc(w)
colorbar
axis off
axis square
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
ph4 = p(18);
f41 = p(19);
f42 = p(20);
f43 = p(21);

g1 = f11.*gampdf(x,f12,f13)./max(gampdf(0:0.01:maxX,f12,f13));
g2 = f21.*gampdf(x,f22,f23)./max(gampdf(0:0.01:maxX,f22,f23));
g3 = f31.*gampdf(x,f32,f33)./max(gampdf(0:0.01:maxX,f32,f33));
g4 = f41.*gampdf(x,f42,f43)./max(gampdf(0:0.01:maxX,f42,f43));

m = g1.*cosd(angle+ph1) + ...
    g2.*cosd(angle.*2+ph2) + ...
    g3.*cosd(angle.*4+ph3) + ...
    g4.*cosd(angle.*8+ph4);

density = (m+1).*(a.*exp(b.*x)+c.*exp(d.*x)+e);

end

