

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

armFilterPointsDegrees = [1.6, 2.1, 4, 6];
armFilterDensityThresh = [9000, 3000, 2000, 1500];


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
    % Filter out unreasonable values from beyond the armFilterPointDegrees
    for ff = 1:length(armFilterPointsDegrees)
        filterIdx = find(supportDeg > armFilterPointsDegrees(ff),1);
        filterRegion = y(:,filterIdx:end);
        filterRegion(filterRegion(:)>armFilterDensityThresh(ff))=nan;
        y(:,filterIdx:end) = filterRegion;
    end
    % Filter out any negative values
    y(y<0)=nan;
    dataMat(:,:,ss)=y;
end

%% Initial fit to mean density profile

% Density, eccentricity, polar angle
Y = nanmean(nanmean(dataMat,3));
X = supportDeg;
P = zeros(size(X));

% objective
maxSupportDeg = supportDeg(find(~isnan(Y), 1, 'last' ));
validIdx = ~isnan(Y);
myObj = @(p) norm(Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p));

% p0 and bounds
asymptoteIdx = find(supportDeg>10,1);
asymptote = nanmean(Y(asymptoteIdx:end));

mBlock0 = [0 0 1.01 1.01];

p0 = [1300, -0.25, 8500, -1.25, asymptote, repmat(mBlock0,1,3)];
lb = [0,-5,0,-5,asymptote, repmat(mBlock0,1,3)];
ub = [1e4,0,1e4,asymptote,1e4, repmat(mBlock0,1,3)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub);


%% Fit polar angle variation

% Density, eccentricity, polar angle
Y = nanmean(dataMat,3);
X = repmat(supportDeg,supportLength,1);
P = repmat(linspace(0,360,supportLength)',1,supportLength);

% objective
validIdx = ~isnan(Y);
myObj = @(p) norm(Y(validIdx) - myModel(X(validIdx),P(validIdx),maxSupportDeg,p));

% p0 and bounds
mBlockLB = [-180 -1 1.01 1.01];
mBlockUB = [180 1 20 20];

p0 = [p(1:5), repmat(mBlock0,1,3)];
lb = [p(1:5), repmat(mBlockLB,1,3)];
ub = [p(1:5), repmat(mBlockUB,1,3)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub);

% generate the model fit
Yfit = nan(supportLength,supportLength);
Yfit(:,:)=myModel(X,P,maxSupportDeg,p);

% plot
figure
surf(X,P,Yfit,'FaceAlpha',0.5,'EdgeColor','none')
hold on
plot3(X(:),P(:),Y(:),'.k')
meridianLabels = {'Nasal','Superior','Temporal','Inferior','Nasal'};
meridianAngles = [0 90 180 270 360];
yticks(meridianAngles);
yticklabels(meridianLabels);
xlabel('Eccentricity [deg]');
zlabel('Density [cones/deg^2]');


figure
for ii = 10:25:110
    semilogy(Y(:,ii),'.');
    hold on
    semilogy(Yfit(:,ii),'-r');
    text(300,Yfit(45,ii),sprintf('%2.2f°',supportDeg(ii)));
end
xticks(meridianAngles);
xticklabels(meridianLabels);
ylabel('log_1_0 density [cones/deg^2]')
foo=1;

figure
for mm=1:4
    subplot(2,2,mm)
    plot(supportDeg,Y(meridianAngles(mm)+1,:),'.k');
    hold on
    plot(supportDeg,Yfit(meridianAngles(mm)+1,:),'-r');  
    xlabel('Eccentricity [deg]');
    ylabel('Density [cones/deg^2]');
    title(meridianLabels{mm});
end

figure
meridianSpec = {'-r','--b','-m','--g'};
for mm=1:4
    plot(supportDeg,Yfit(meridianAngles(mm)+1,:),meridianSpec{mm},'LineWidth',2);
    hold on
end
xlabel('Eccentricity [deg]');
ylabel('Density [cones/deg^2]');
legend(meridianLabels(1:4));





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

