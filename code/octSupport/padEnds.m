function paddedLoc = padEnds(inLoc,padSize)

[MaxVal, maxI] = max(inLoc(:,1));
[MinVal, minI] = min(inLoc(:,1));

paddedLoc = [inLoc; [(MaxVal+1):(MaxVal+padSize)]' ones(padSize,1)*inLoc(maxI,2); ...
    [(MinVal-padSize):(MinVal-1)]' ones(padSize,1)*inLoc(minI,2)];

