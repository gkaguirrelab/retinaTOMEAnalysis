function saveSingleVerticalOCT(dataDir,saveDir)
%Purpose: Montages the extended horizontal OCT into a single aligned image
%Inputs:
%dataDir - Input directory with all OCT data
%left/center/right images in the montage
%saveDir - output directory to save the montaged images

FullSave = 0;%flat for saving additional information about the alignments
allSubs = dir(fullfile(dataDir,'1*'));%find all subject directories
eyeSide = {'OS', 'OD'};
counter = 0;
%Process for all subjects and all eyes
for s = 1:length(allSubs)
    for eye = 1:2
        %Display current subject info
        subID=allSubs(s).name
        eyeSide{eye}
        %set current subject input and save directories
        inDirSub=fullfile(dataDir,subID,'HeidelbergSpectralisOCT',eyeSide{eye});
        saveDirSub = fullfile(saveDir,subID,eyeSide{eye});
        mkdir(saveDirSub);%make it if it's not made already
        
        %find all vol files in directory
        allVols = dir(fullfile(inDirSub,'*.vol'));
        
        %search data for vertical scans
        for i = 1:length(allVols)
            PATH = fullfile(inDirSub,allVols(i).name);
            [HEADER, BSCANHEADER, SLO, Bscans] = openVolFast(PATH,'nodisp');
            
            %Check for vertical and and single images
            if(HEADER.NumBScans == 1  && sum(abs(BSCANHEADER.StartX - BSCANHEADER.EndX)) < 1)
                counter = counter+1;
                Bscans = Bscans.^.25;
                imOut = Bscans(:,:,1);
                res = [HEADER.ScaleZ HEADER.ScaleX .001];
                
                %export as nii for delineation
                origin = [0 0 0];
                %Reorient matlab img to nii img and save
                img = permute(imOut,[2 1 3]);
                img = rot90(img,2);
                niires = [res(2) res(1) res(3)];
                nii = make_nii(img, niires, origin);
                save_nii(nii, fullfile(saveDirSub,[subID '_' eyeSide{eye} '_SingleVerticalOCT.nii']))
            end
        end
        
     end
end
counter