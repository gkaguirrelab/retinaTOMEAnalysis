function filteredIm = findOrderedRegions(im,precision)

% Specify the precision if not defined
if nargin < 2
    precision = -3;
end

% Copy over im to roundIm;
roundIm = im;

% Set nan values in im to -1
roundIm(isnan(roundIm))=-1;

% Round the image values to integers with the specified precision. A value
% of -2, for example, rounds to the nearest value of 100.
roundIm = round(roundIm,precision);

% Get the ordered list of pixel values
pixValList = sort(unique(roundIm(:)),'descend');

% Set the initial region to be the maximum pixel value
lastOrderedRegions = {find(roundIm == pixValList(1))};
lastBorderRegions = lastOrderedRegions;
lastMaxBorderVals = pixValList(1);

% Set the start point of the search to the second value
hIdx = 2;

% Search down through H
stillSearching = true;
while stillSearching
    
    % Check if we have accounted for all the pixels
    if hIdx > length(pixValList)
        stillSearching = false;
        hIdx = hIdx - 1;
        continue
    end
    
    % Get the current threshold
    H = pixValList(hIdx);
    
    % Binarize the image at this threshold
    binIm = roundIm;
    binIm(roundIm<H)=0;
    binIm(roundIm>=H)=1;
    
    % Find the connected regions, and sort the regions by max value
    CC = bwconncomp(binIm);
    regionMax = cellfun(@(x) max(roundIm(x)),CC.PixelIdxList);
    [sortList,sortIdx]=sort(regionMax,'descend');
    sortIdx = sortIdx(~isnan(sortList));
    orderedRegions = CC.PixelIdxList(sortIdx);
    
    % Get the number of new regions
    nNewRegions = length(orderedRegions) - length(lastOrderedRegions);
    
    % If the number of regions has decreased, then we are done searching
    if nNewRegions < 0
        stillSearching = false;
        hIdx = hIdx - 1;
        continue
    end
    
    % Compare the current set of orderedRegions to the last set
    % of regions
    newBorderRegions = {};
    for rr = 1:length(lastOrderedRegions)
        
        % Find the difference between the region now and the region at the
        % prior threshold
        borderRegion = setdiff(orderedRegions{rr},lastOrderedRegions{rr});
        
        % If the region didn't change in size, retain the prior border.
        % Otherwise, retain just the area that the border has grown
        if isempty(borderRegion)
            newBorderRegions{rr} = lastBorderRegions{rr};
        else
            newBorderRegions{rr} = setdiff(orderedRegions{rr},lastOrderedRegions{rr});
        end
    end
    
    % If there are new regions, add these entire regions as new
    % borders
    if nNewRegions > 0
        newBorderRegions(end+1:end+nNewRegions) = orderedRegions(end+1-nNewRegions:end);
    end
    
    % Obtain the max value in the new border regions
    maxBorderVals = cellfun(@(x) max(roundIm(x)),newBorderRegions);
    
    % Check if any of these values are larger than in the prior border
    % regions. If so, then we are done our search. If not, move to the next
    % pixel value.
    if any(maxBorderVals(1:length(lastMaxBorderVals)) > lastMaxBorderVals)
        stillSearching = false;
        hIdx = hIdx - 1;
        continue
    else
        lastOrderedRegions = orderedRegions;
        lastBorderRegions = newBorderRegions;
        lastMaxBorderVals = maxBorderVals;
        hIdx = hIdx + 1;
    end
    
end

% Generate the filtered image

% The final threshold
H = pixValList(hIdx);

% Binarize the image at this threshold
binIm = roundIm;
binIm(roundIm<H)=0;
binIm(roundIm>=H)=1;

% Find the connected regions
CC = bwconncomp(binIm);

% Generate the image
filteredIm = nan(size(roundIm));

% Move the connected objects from the original im to the filtered image
for rr = 1:CC.NumObjects
    filteredIm(CC.PixelIdxList{rr}) = im(CC.PixelIdxList{rr});
end


end
