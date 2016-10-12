function [sapPts] = readMAIA(xmlFile)
%
% [sapPts] = readMAIA(xmlFile)
% This function reads in the xml file associated with the MAIA perimetry
% data and turns it into a matrix.
%
% Inputs:
% xlmFile = the projection.xlm file located in the MAIA directory (as of now they are in zipped folder)
% 
% Outputs:
% sapPts = an mx4 matrix of perimitry data. 
%          Column order: [id, ray, angle_deg, final_intensity]
%                 Units: [#, deg, deg, dB] 
% MAB OCT 2016

doc = xml2struct(xmlFile);

for i = 1:length(doc)
   if strcmp(doc(i).Name,'MitProjection')
       idx1 = i;
   end
end

for j = 1:max(size(doc(idx1).Children))
    if strcmp(doc(idx1).Children(j).Name,'Stimuli')
       idx2 = j;
    end
end

count = 0;
for k = 1:max(size(doc(idx1).Children(idx2).Children))
    if strcmp(doc(idx1).Children(idx2).Children(k).Name,'Stimulus')
        count = count+1;
        for m = 1:length(doc(idx1).Children(idx2).Children(k).Attributes)
            if strcmp(doc(idx1).Children(idx2).Children(k).Attributes(m).Name,'id')
                sapPts(count,1) = str2num(doc(idx1).Children(idx2).Children(k).Attributes(m).Value);
            elseif strcmp(doc(idx1).Children(idx2).Children(k).Attributes(m).Name,'ray')
                sapPts(count,2) = str2num(doc(idx1).Children(idx2).Children(k).Attributes(m).Value);
            elseif strcmp(doc(idx1).Children(idx2).Children(k).Attributes(m).Name,'angle_deg')
                sapPts(count,3) = str2num(doc(idx1).Children(idx2).Children(k).Attributes(m).Value);
            elseif strcmp(doc(idx1).Children(idx2).Children(k).Attributes(m).Name,'final_intensity')
                sapPts(count,4) = str2num(doc(idx1).Children(idx2).Children(k).Attributes(m).Value);
            end
        end
    end
end



end
