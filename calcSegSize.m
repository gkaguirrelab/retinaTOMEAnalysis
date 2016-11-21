function areaPerSeg= calcSegSize(radDeg,smpPerDeg,sectorAngle)

binWidthEcc =1/smpPerDeg;
r_inner = 0:binWidthEcc:radDeg-binWidthEcc;
r_outer = binWidthEcc:binWidthEcc:radDeg;

areaPerSeg =  ((pi.*r_outer.^2)-(pi.*r_inner.^2))./(360/sectorAngle);
areaPerSeg = [0 areaPerSeg];
end