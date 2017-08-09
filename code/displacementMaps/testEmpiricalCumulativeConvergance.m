% #############################################################
% The goal is to produce a plot of the cumulative retinal
% ganglion cell and receptive field counts as a function of
% eccentricity (mm).
%
% MAB 2016
%
% Notes:
% THIS CODE SETS THE REFERENCE FRAME TO THE RETINA
%
% #############################################################

% housekeeping
clc
clear vars
close all

% Define the meridia for which we will make the calculation
meridia = [0 90 180 270];
meridiaNames = {'nasal' 'superior' 'temporal' 'inferior'};
displacementPointDeg = [17 17 11 17];
midgetFracAdjust = [0.8230 1.0375 0.8453 1.1791];

% Loop over meridia
for ii = 1:length(meridia)
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ midget_rgcDensity_degSq, supportPosDeg ] = getCurcioMidgetRGCDensityByEccen( meridia(ii) );
    
    % Plot the raw RGC density data
    plot(supportPosDeg,midget_rgcDensity_degSq,'.r')
    xlabel('eccentricity [deg]');
    ylabel('density [cells / deg2]');
    title('midget RGC density derived from Curcio and Allen 1990 (w spline fit)');
    hold on
    
    % Adjust the midget RGC denisty by the fudge factor
    midget_rgcDensity_degSq = midget_rgcDensity_degSq * midgetFracAdjust(ii);

    % Fit a spline to the RGC density data
    splineFunctionSuperior = fit(supportPosDeg,midget_rgcDensity_degSq,'smoothingspline', 'Exclude',find(isnan(midget_rgcDensity_degSq)),'SmoothingParam', 1);
    
    % Add the spline fit to the plot
    plot(0:0.1:max(supportPosDeg), splineFunctionSuperior(0:0.1:max(supportPosDeg)),'-k');
    hold off
    
    % Switch to equal eccentricity sampling between 0 and 20 degrees
    supportPosDeg = 0:0.01:30;
    
    % Now get the midget RF density
    midgetRFDensity_degSq = getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, 90);
        
    % Plot the midget RF density
    figure
    plot(supportPosDeg,midgetRFDensity_degSq,'-b')
    xlabel('eccentricity [deg]');
    ylabel('density [mRF / deg2]');
    title('midget RF density from Watson 2014 eq 8');

    % Calculate the counts per growing ring
    
    % Watson 2*pi*r correction
    ringArea = [0,diff(supportPosDeg.^2 * pi)];
    
    countPerRingRF = cumsum(midgetRFDensity_degSq.*ringArea);
    countPerRingRGC = cumsum(splineFunctionSuperior(supportPosDeg).*ringArea');
    
    % Plot the RGC and RF density data, within growing areas
    figure
    plot(supportPosDeg,countPerRingRGC,'-k')
    xlabel('eccentricity [deg]');
    ylabel('counts [cells per sector]');
    hold on
    plot(supportPosDeg,countPerRingRF,'-b')
    title([meridiaNames{ii} ' - counts per ring, mRGC black, mRF blue']);
    hold off
    
    % Calculate the RF / RGC ratio within the 20-30 degree range
%     idx = find((supportPosDeg > 20) .* (supportPosDeg < 30));
%     mean(countPerRingRF(idx) ./ countPerRingRGC(idx)')

    % Calculate the RF / RGC ratio at the theoretical displacement point
    idx = find(supportPosDeg == displacementPointDeg(ii));
    mean(countPerRingRF(idx) ./ countPerRingRGC(idx)')
        
    pause
    close all
    
end % loop over meridia

