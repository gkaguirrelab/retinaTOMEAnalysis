
function countsPerRing = calcCumulative(regularSupportPosDeg, densityFunction)

ringArea = [0,diff(regularSupportPosDeg.^2 * pi)];
countsPerRing = cumsum(densityFunction.*ringArea);

end