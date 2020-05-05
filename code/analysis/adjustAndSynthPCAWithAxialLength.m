nDimsToUse = 6;

corrCoeffAxialLength = nan(1,nDimsToUse);
fprintf('Correlation axial length with PCA coefficients:\n');
for ii = 1:nDimsToUse
    corrCoeffAxialLength(ii) = corr(axialLengths',coeff(:,ii));
    txt = sprintf('\tPC%d: %2.2f \n',ii,corrCoeffAxialLength(ii));
    fprintf(txt);
end

% Regress the axial length upon the coefficients
X = [coeff(:,1:nDimsToUse)'; ones(1,50)]';
b = X\axialLengths';

% % Show regression fit of PCA and axial length
figure
plot(axialLengths,X*b,'*r');
axis square
refline(1,0);


adjustedCoeff = zeros(size(coeff));
ALmax= max(axialLengths);
ALmin= min(axialLengths);
ALRange = linspace(ALmin,ALmax,50);
ALRange(16) = 23.58;%We're forcing this one to be the mean emmetropic value (was 23.56 before).
synCoeff = zeros(50,nDimsToUse);
for p=1:nDimsToUse
    
    currCoeff = coeff(:,p);
    %find regression line between PC coeff and  axial length
    pp = polyfit(axialLengths,currCoeff',1);
    
    %find PC coeff value for emmetrope
    %Adjust each score by fraction: coeff_emm / coeff_regress
    AL_emmetrope_mm=23.58;
    adjustProportion_allSub = polyval(pp,AL_emmetrope_mm) - polyval(pp,axialLengths);
    adjustedCoeff(:,p) = currCoeff + adjustProportion_allSub';
    
    %Save synthesized representation of AL
    synCoeff(:,p) = polyval(pp,ALRange);
end