function midgetFraction = calcWatsonMidgetFractionByEccen(supportPosDeg,f0,rm)

% The equation is taken from Watson JoV 2014 (eq 7, plotted in fig 8), 
% which Watson took from Drasdo et al 2007.
midgetFraction = f0.*(1+(supportPosDeg./rm)).^-1;

end