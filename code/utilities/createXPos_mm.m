% This function creates a matrix of XPos_mm, where each subject has their
% support function in mm, which is going to vary over subjects principally
% due to differences in the axial length of the eye


% Where is DropBox?
dropboxBaseDir=fullfile(getpref('retinaTOMEAnalysis','dropboxBaseDir'));

% Load the mmPerDeg polynomial conversion functions
mmPerDegFileName = fullfile(dropboxBaseDir,'AOSO_analysis','mmPerDegMaps','mmPerDegPolyFit.mat');
load(mmPerDegFileName)

% Loop through the horizontal and vertical meridians
for mm = 1:2
    
    % Load the appropriate data
    switch mm
        case 1
            orientation = 'horiz';
            GCIPThickFileName = fullfile(dropboxBaseDir,'AOSO_analysis','OCTExplorerExtendedHorizontalData','GCIP_thicknessesByDeg.mat');
            load(GCIPThickFileName)
        case 2
            orientation = 'vert';
            GCIPThickFileName = fullfile(dropboxBaseDir,'AOSO_analysis','OCTSingleVerticalData','GCIP_thicknessesByDeg.mat');
            load(GCIPThickFileName)
    end
    
    % Loop through the subjects
    for ii=1:size(GCthicknessValuesAtXPos_um,1)
        
        pp = mmPerDegPolyFit{ii};
        switch orientation
            case 'horiz'
                mmSqPerDegSq = pp([-XPos_Degs;zeros(size(XPos_Degs))]').^2;
            case 'vert'
                mmSqPerDegSq = pp([zeros(size(XPos_Degs));-XPos_Degs]').^2;
        end
        
        % Create the XPos_mm for each subject
        eccenPos = @(hh) pp(hh, 0);
        for xx = 1:length(XPos_Degs)
            if XPos_Degs(xx)==0
                continue
            end
            if XPos_Degs(xx)>0
                XPos_mm(xx,ii,mm) = integral(eccenPos,0,XPos_Degs(xx));
            else
                XPos_mm(xx,ii,mm) = -integral(eccenPos,XPos_Degs(xx),0);
            end
        end
    end
    
end