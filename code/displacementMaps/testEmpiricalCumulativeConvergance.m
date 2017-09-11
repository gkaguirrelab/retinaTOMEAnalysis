% #############################################################
% The goal is to calculate RGC displacement as a function of eccentricity
% for the four cardinal merdians, based upon the assumption
% that the cumulative counts of mRGCs and mRFs are lined across the retina,
% following the approach of Watson 2014 and Turpin et al. 2015.
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

% define some hard-coded params
sampleResolutionDegrees = 0.01;

% Define the meridians for which we will make the calculation
meridians = [0 90 180 270];
meridianNames = {'Nasal' 'Superior' 'Temporal' 'Inferior'};
displacementPointDeg = [17 17 11 17];
midgetFracAdjust = [0.8230 1.0375 0.8453 1.1791];

% Loop over meridians
for ii = 1:length(meridians)
    
    % Load the RGC Density Data from Curcio and Allen 1990:
    [ midget_rgcDensity_degSq, supportPosDeg ] = getCurcioMidgetRGCDensityByEccen( meridians(ii) );

    % Adjust the midget RGC denisty by the fudge factor
    midget_rgcDensity_degSq = midget_rgcDensity_degSq * midgetFracAdjust(ii);
    
    % Fit a spline to the RGC density data
    splineFunctionSuperior = fit(supportPosDeg,midget_rgcDensity_degSq,'smoothingspline', 'Exclude',find(isnan(midget_rgcDensity_degSq)),'SmoothingParam', 1);
    
    % Plot the raw and fitted RGC density data
%     plot(supportPosDeg,midget_rgcDensity_degSq,'.r')
%     xlabel('eccentricity [deg]');
%     ylabel('density [cells / deg2]');
%     title('midget RGC density derived from Curcio and Allen 1990 (w spline fit)');
%     hold on
%     plot(0:0.1:max(supportPosDeg), splineFunctionSuperior(0:0.1:max(supportPosDeg)),'-k');
%     hold off
    
    % Switch to equal eccentricity sampling between 0 and 20 degrees
    supportPosDeg = 0:sampleResolutionDegrees:30;
    
    % Now get the midget RF density
    midgetRFDensity_degSq = getWatsonMidgetReceptiveFieldDensityByEccen(supportPosDeg, 90);
        
    % Plot the midget RF density
%     figure
%     plot(supportPosDeg,midgetRFDensity_degSq,'-b')
%     xlabel('eccentricity [deg]');
%     ylabel('density [mRF / deg2]');
%     title('midget RF density from Watson 2014 eq 8');

    % Calculate the counts per Watson integration ring
    ringArea = [0,diff(supportPosDeg.^2 * pi)];    
    countPerRingRF = cumsum(midgetRFDensity_degSq.*ringArea);
    countPerRingRGC = cumsum(splineFunctionSuperior(supportPosDeg).*ringArea');
    
    % Plot the RGC and RF density data, within growing areas
    figure
    subplot(2,1,1);
    plot(supportPosDeg,countPerRingRGC,'-k')
    xlabel('eccentricity [deg]');
    ylabel('counts [cells per sector]');
    hold on
    plot(supportPosDeg,countPerRingRF,'-b')
    title([meridianNames{ii} ' - counts per ring, mRGC black, mRF blue']);
    hold off
    
    % Calculate the RF / RGC ratio within the 20-30 degree range
%     idx = find((supportPosDeg > 20) .* (supportPosDeg < 30));
%     mean(countPerRingRF(idx) ./ countPerRingRGC(idx)')

    % Calculate the RF / RGC ratio at the theoretical displacement point
%     idx = find(supportPosDeg == displacementPointDeg(ii));
%     mean(countPerRingRF(idx) ./ countPerRingRGC(idx)')

    % Measure the displacement (in degrees)
    displaceInSamples=arrayfun(@(x) find(countPerRingRGC>x,1), countPerRingRF(2:end),'UniformOutput',false);
    displaceInSamples=cell2mat(displaceInSamples(cellfun(@(x) ~isempty(x), displaceInSamples)));
    displaceInDeg = (displaceInSamples - (1:1:length(displaceInSamples))) * sampleResolutionDegrees;
    displaceInDeg(find(displaceInDeg <= 0,1):end)=0;
    
    % Plot the displacement
    subplot(2,1,2);
    plot(supportPosDeg(1:length(displaceInDeg)),displaceInDeg,'-r')
    xlabel('eccentricity [deg]');
    ylabel('radial displacement of RGCs [deg]');
    title([meridianNames{ii} ' - RGC displacement']);
        
end % loop over meridia

