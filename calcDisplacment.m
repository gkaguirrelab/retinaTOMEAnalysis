function disMap = calcDisplacment(Df,dRGC,radDeg,smpPerDeg)

% Add minima finder to center data.

numXsmps = size(Df,2);
numYsmps = size(Df,1);

Xzero = round(numXsmps/2);
Yzero = round(numYsmps/2);

Xsmps = linspace(Xzero,Xzero+smpPerDeg*radDeg,smpPerDeg*radDeg+1);

%calculate the cumulative denstiy count 

for i = 1:length(Xsmps)
    
    if i == 1
        dfCount(1,i) =Df(Yzero,Xsmps(i))
        rgcCount(1,i) =dRGC(Yzero,Xsmps(i))
    else
        dfCount(1,i) = dfCount(1,i-1) + Df(Yzero,Xsmps(i));
        drgcCount(1,i) = drgcCount(1,i-1) + dRGC(Yzero,Xsmps(i));
    end
end


end


