subjects = [11015,11018,11028,11031,11043,11049,11050,11051,11052,11053,11055,11056,11057,11058,11060,11061,11062,11064,11065,11067,11068,11069,11070,11071,11072,11073,11074,11075,11076,11077,11078,11080,11081,11082,11083,11084,11086,11087,11088,11089,11091,11092,11093,11094,11095,11096,11097,11098,11099,11100];
axialLength = [27.52,27.07,24.03,21.89,25.15,24.49,24.9,26.05,22.76,23.89,27.49,25.1,22.18,22.96,23.91,24.92,24.79,22.69,23.45,22.9,24.83,23.54,25.35,22.64,26.54,25.47,25.29,24.18,25.35,24.66,25.35,22.44,23.45,22.88,23.67,23.38,23.73,24.01,24.06,24.5,25.7,26.57,23.56,23.39,22.17,26.24,23.08,24.23,22.11,24.85];
sphericalAmetropia = [-7.5,-4.75,0.25,3.5,-1.75,-3.25,-3.75,-8.5,-0.5,-5.25,-10.25,-3.25,0.5,-0.25,-0.75,-2,-2,-1.25,-0.5,0.75,-1,-0.5,-0.5,0,-6,-2,-5.25,-0.75,-5,-4.25,-1.5,0.5,-0.75,0.25,0.5,0,-0.25,-1.75,-0.25,0.25,-5.25,-5.25,-1,-0.5,-1,-6.25,-0.25,-1.5,0.25,-5];
k1 = [41.56,nan,43.38,43.38,41.77,43.32, 44.2,nan,45.30, 45.7,41.62,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,42.08,nan,nan,41.36,43.49,44.35,45.42,42.61,43.77,42.67,45.49,42.03,41.21,41.87,42.61,45.61,43.66,46.87,43.38,42.45,42.72,45.92,44.00];
k2 = [42.35,nan,43.55,43.95,43.55, 44.4,44.94,nan,46.23, 46.6,43.66,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,42.83,nan,nan,nan,44.64,45.36,45.42,43.72,44.70,43.05,45.86,42.51, 42.2,42.72,43.89,45.86,44.64,47.34,43.89,43.55,45.06,47.14,45.30];
angle = [5,nan,153,166,25,23,156,nan,8,0,174,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,16,nan,nan,25,3,4,42,173,175,18,166,1,178,26,28,48,178,163,165,8,172,17,5];

saveDir = '~/Dropbox (Aguirre-Brainard Lab)/AOSO_analysis';

resultSet = {};

parfor (ii = 1:2, 2) %length(subjects)

    cc=[k1(ii),k2(ii),angle(ii)];
    if any(isnan(cc))
        cc = [];
    end
    SR = sphericalAmetropia(ii);
    if isnan(SR)
        SR=[];
    end
    eye = modelEyeParameters('axialLength',axialLength(ii),'sphericalAmetropia',SR,'measuredCornealCurvature',cc);
    
    % Create a matrix of 
    S = eye.retina.S;
    alpha = [5.45 2.5 0];
    horizVals = -15:1:15;
    vertVals = -15:1:15;

    mmSqPerDegSq = nan(length(horizVals),length(vertVals));

    for jj = 1:length(horizVals)
        for kk = 1:length(vertVals)
            degField = [horizVals(jj) vertVals(kk) 0] + alpha;
            [~,X0] = findRetinaFieldPoint( eye, degField);
            degField = degField + [0.0707 0.0707 0];
            [~,X1] = findRetinaFieldPoint( eye, degField);
            distance = quadric.panouGeodesicDistance(S,[],[],X0,X1);
            mmSqPerDegSq(jj,kk) = (distance*10)^2;
        end
    end
    fprintf(['Done subject ' num2str(subjects(ii))]);
    resultSet(ii) = {mmSqPerDegSq};
end

for ii = 1:length(subjects)
    outfile = fullfile(saveDir,'mmSqPerDegSqMaps',[num2str(subjects(ii)) '_mmSqPerDegSqMap.mat']);
    mmSqPerDegSq = resultSet{ii};
    save(outfile,'mmSqPerDegSq');
end
