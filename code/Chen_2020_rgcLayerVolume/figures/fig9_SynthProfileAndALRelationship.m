function fig9_SynthProfileAndALRelationship(GCVolPCACoeff,axialLengths,XPos_Degs,comboTable,scoreExpandedSmoothed,nDimsToUse,saveDir,orientation)

% Plot the synthesized reconstructions by axial length
minAL = min(comboTable.Axial_Length_average);
maxAL = max(comboTable.Axial_Length_average);
emAL = 23.58;
ALs = [minAL emAL maxAL];

for ii=1:length(ALs)
    for dd = 1:nDimsToUse
        
        % find regression line between PC coeff and  axial length
        pp = polyfit(axialLengths,GCVolPCACoeff(:,dd),1);
        coeff(dd) = polyval(pp,ALs(ii));
    end
    
    synthProfileVol(:,ii) = scoreExpandedSmoothed(:,1:nDimsToUse)*coeff';
    
end


% Plot the GC volume functions
str = sprintf('Synthesized GC vol profiles for AL = %2.2f, %2.2f, %2.2f',ALs);
h=profilePlot(XPos_Degs, synthProfileVol, [], 'Eccentricity [deg visual angle]','GC Tissue Volume [mm^3 / deg^2]',[],1);
title(str)
ylim([0 8e-3]);
xlim([-25 25]);
setTightFig
saveas(h,fullfile(saveDir,'fig9','a.pdf'));

% Make a mmPerDeg function for each of these model eyes, and produce the
% thickness profiles, and thickness profiles by mm
XPos_mm = zeros(length(XPos_Degs),length(ALs));

% Create the down-sampled XPos_Degs that we will measure
XPos_Support = linspace(1,length(XPos_Degs),100);
XPos_DegSub = interp1(1:length(XPos_Degs),XPos_Degs,XPos_Support);

% Loop over model eyes
for ss = 1:length(ALs)
    eye = modelEyeParameters('axialLength',ALs(ss));
    foveaCoord = eye.landmarks.fovea.coords';
    fieldAngularPosition = eye.landmarks.fovea.degField;
    rayOriginDistance = 1500;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    distanceReferenceCoord = calcPrincipalPoint(eye);
    
    mmPerDeg = nan(size(XPos_DegSub));
    retinalDistance = nan(size(XPos_DegSub));
    for ii = 1:length(XPos_DegSub)
        switch orientation
            case 'horiz'
                deltaAngles = [1 0]; % Measure horizontally separated points
                thisPosition = fieldAngularPosition + [ -XPos_DegSub(ii) 0 ];
            case 'vert'
                deltaAngles = [1 0]; % Measure vertically separated points
                thisPosition = fieldAngularPosition + [ 0 -XPos_DegSub(ii) ];
        end
        rayPath = calcNodalRayFromField(eye,thisPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);
        retinalDistance(ii) = sign(XPos_DegSub(ii))*quadric.geodesic(eye.retina.S,[foveaCoord,rayPath(:,end)]);
        mmPerDeg(ii) = calcMmRetinaPerDeg(eye,thisPosition,deltaAngles,rayOriginDistance,angleReferenceCoord);

    end
    mmSqPerDegSq(:,ss) = interp1(XPos_DegSub,mmPerDeg,XPos_Degs).^2;
    XPos_mm(:,ss) = interp1(XPos_DegSub,retinalDistance,XPos_Degs);
    synthProfileThick(:,ss) = synthProfileVol(:,ss)./mmSqPerDegSq(:,ss);
end


% Plot the GC thick functions with support in degrees
str = sprintf('Synthesized GC thickness profiles for AL = %2.2f, %2.2f, %2.2f',ALs);
h=profilePlot(XPos_Degs, synthProfileThick, [], 'Eccentricity [deg visual angle]','GC Tissue Thickness [mm]',[],1);
ylim([0 0.075]);
xlim([-25 25]);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig9','b.pdf'));

% Save a data file that contains the profile of the GC thickness in mm as a
% function of eccentricity in degrees
switch orientation
    case 'horiz'
        thicknessMmHoriz = squeeze(synthProfileThick(:,2));        
        save(fullfile(saveDir,'emmetropeThickProfile_horiz.mat'),'XPos_Degs','thicknessMmHoriz');
    case 'vert'
        thicknessMmVert = squeeze(synthProfileThick(:,2));        
        save(fullfile(saveDir,'emmetropeThickProfile_vert.mat'),'XPos_Degs','thicknessMmVert');
end

% Plot the GC thick functions with support in mm
h = figure;
str = sprintf('Synthesized GC thickness profiles for AL = %2.2f, %2.2f, %2.2f',ALs);
for ii = 1:length(ALs)
    plot(XPos_mm(:,ii), synthProfileThick(:,ii), '-');
    hold on
end
xlabel('Eccentricity [mm]');ylabel('GC Tissue Thickness [mm]');
ylim([0 0.075]);
xlim([-8 8]);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig9','c.pdf'));

end


