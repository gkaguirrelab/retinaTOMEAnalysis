function fig9_SynthProfileAndALRelationship(GCVolPCACoeff,axialLengths,XPos_Degs,comboTable,scoreExpandedSmoothed,nDimsToUse,saveDir)
% Plot the synthesized reconstructions by axial length
[minAL,hyperopeIdx] = min(comboTable.Axial_Length_average);
[maxAL,myopeIdx] = max(comboTable.Axial_Length_average);
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
setTightFig
saveas(h,fullfile(saveDir,'fig9','a.png'));

% Make a mmPerDegMap for each of these model eyes, and produce the
% thickness profiles, and thickness profiles by mm
for ii=1:length(ALs)
    mmPerDegPolyFit{ii} = magMap(ALs(ii));
    mmSqPerDegSq = mmPerDegPolyFit{ii}([zeros(size(XPos_Degs));-XPos_Degs]').^2;
    synthProfileThick(:,ii) = synthProfileVol(:,ii)./mmSqPerDegSq;

end

% Plot the GC thick functions
str = sprintf('Synthesized GC thickness profiles for AL = %2.2f, %2.2f, %2.2f',ALs);
h=profilePlot(XPos_Degs, synthProfileThick(:,idx), [], 'Eccentricity [deg visual angle]','GC Tissue Thickness [mm]',[],1);
title(str)
setTightFig
saveas(h,fullfile(saveDir,'fig9','b.png'));




end




function mmPerDegPolyFit = magMap(AL)

% Create the model eye
eye = modelEyeParameters('axialLength',AL,'calcLandmarkFovea',true);

% Define the visual field domain over which we will make the measure,
% in degrees relative to the fovea
horizVals = -30:15:30;
vertVals = -30:15:30;

% Define an empty matrix to hold the results
mmPerDeg = nan(length(horizVals),length(vertVals));

% Define the delta deg
deltaDegEuclidean = 1e-3;
deltaAngles = [sqrt(deltaDegEuclidean/2) sqrt(deltaDegEuclidean/2)];

% Loop over horizontal and vertical field positions
for jj = 1:length(horizVals)
    for kk = 1:length(vertVals)
        % The position in the field relative to the optical axis of the
        % eye
        
        degField = [horizVals(jj) vertVals(kk)] + eye.landmarks.fovea.degField(1:2);
        
        % Obtain the retinal points that are delta degrees on either
        % side of the specified degree field position
        [~,X0] = calcRetinaFieldPoint( eye, degField - deltaAngles./2);
        [~,X1] = calcRetinaFieldPoint( eye, degField + deltaAngles./2);
        
        % The difference between X0 and X1 is used to calculate the mm
        % of retina per degree of visual field for this location
        mmPerDeg(jj,kk) = norm(X0-X1) / norm(deltaAngles);
        
    end
    
end

% Fit a polynomial surface to the measure
[X,Y]=meshgrid(horizVals,vertVals);
mmPerDegPolyFit = fit([X(:),Y(:)],mmPerDeg(:),'poly33');


end