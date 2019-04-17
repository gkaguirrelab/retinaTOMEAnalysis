% Calc RGC tissue volume

% Best corrected peripheral acuity varies minimally with spherical
% refractive error, when this acuity is expressed in cycles / deg and
% is not influenced by the minimizing effect of spectacle lenses. The most
% direct measurement of this is with laser interference fringes that are
% presented directly to the retina.
%
% Support for this claim is found in:
%
%   Coletta, Nancy J., and Tonya Watson. "Effect of myopia on visual acuity
%   measured with laser interference fringes." Vision research 46.5 (2006):
%   636-651.
%
%   Shown in Figure 3, there is a minimal, non-significant dependence upon
%   foveal and periperal acuity as a function of spherical refractive
%   error.
%
%   Atchison, David A., Katrina L. Schmid, and Nicola Pritchard. "Neural
%   and optical limits to visual performance in myopia." Vision Research
%   46.21 (2006): 3707-3722.
%
%   Table 1 shows a small, non-significant change in peripheral acuity with
%   spherical refractive error as measured with laser interference in 121
%   subjects.
%
% These data are accounted for by a model of "retinal stretching" in which
% there is common endowment of retinal ganglion cells across people, but
% that these cells are distributed uniformly across an eye that has
% undergone expansion of the posterior pole in myopia.
%
% These findings suggest that, in a population of young, healthy subjects,
% a measurement of the retinal ganglion cell density per degree square of
% the visual field should not differ by spherical ametropia. Measurements
% of RGC+IPL layer thickness by OCT, however, are found to be strongly
% negatively correlated with axial length. This may be explained by the
% retinal stretching model.
%
% To account for this, we can convert measurements of retinal thickness
% into measurements of tissue volume. To do so, we obtain the mm of retina
% per degree of visual angle based upon the axial length and spherical
% refractive error of the person.
%
% The conversion assumes that there is some fixed amount of the RGC+IPL
% layer that is of static thickness, and does not vary by axial length. I
% find that subtracting 45 microns of thickness from the mean RGC+IPL
% thickness, and then multiplying by retinal surface area (mm2 per deg2)
% removes the effect of axial length from the "corrected" RGC tissue
% volume. In our data, the emmetropic eye has a mean RGC+IPL thickness of
% 61 microns. This implies that 25% of the thickness of the RGC+IPL layer
% is RGC cell bodies, which accords well with the finding that RGC
% thickness is well modeled by having 50% of the volume be RGCs, plus
% having an IPL layer that is 40% of the total thickness.


%{
    % Relation between axial length and mm per deg at the retinal apex
    deltaAngles=[sqrt(1/2)/2 sqrt(1/2)/2 0];
    length = [];
    mmPerDeg = [];
    for SR = -7:1:3
        eye = modelEyeParameters('sphericalAmetropia',SR);
        [~,X0] = findRetinaFieldPoint( eye, -deltaAngles);
        [~,X1] = findRetinaFieldPoint( eye, deltaAngles);
        length(end+1) = eye.meta.axialLength;
        mmPerDeg(end+1) = sqrt(sum((X0-X1).^2)) ./ sqrt(sum((deltaAngles.*2).^2));
    end
    figure
    plot(length,mmPerDeg,'xr');
    xlabel('axial length [mm]');
    ylabel('mm retina per deg visual angle');
    vals = polyfit(length,mmPerDeg,1);
    fprintf('Retinal mm per deg visual field at the viterous chamber apex = (%2.3f * axialLength) %2.3f \n',vals(1),vals(2));
%}

% This is the mmPerDeg at the ellipsoidal pole of the vitreous chamber
mmPerDeg = @(axialLength) (0.0165.*axialLength)-0.1070;

correctedThickness = (thick-45).*mmPerDeg(length).^2;

