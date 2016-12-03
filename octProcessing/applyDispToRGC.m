function correctedRGC = applyDispToRGC(rgcThickness,dispMap)

xCent = round(size(rgcThickness,2)/2);
yCent = round(size(rgcThickness,1)/2);

%%for a point in rcgThickness, look up in dispMap the displacement.
for y = 1:size(rgcThickness,1)
    for x = 1:size(rgcThickness,2)
        
%% find the x y coord of the new point that is disp mm from starting point towards image center

%% iterativly nansum

%% outputs
    end
end


end