%This is a script that goes through all the .vol files we have and view each subject
%Used for manually determining which vol file corresponds to which type of
%scan
%dropboxBaseDir = getpref('retinaTOMEAnalysis','dropboxBaseDir');
inDir=fullfile(dropboxBaseDir,'AOSO_data','connectomeRetinaData');

%change this number to start at different subject numbers
StartingSubjectNumber = 1;

%find all subject directories
allSubs = dir(fullfile(inDir,'1*'));
eyeSide = {'OS', 'OD'};
%process for each subject and each eye
for s = StartingSubjectNumber:length(allSubs)
    for eye = 1:2
        %Display subject and eye
        subID=allSubs(s).name
        eyeSide{eye}
        %current working directory 
        inDirSub=fullfile(inDir,subID,'HeidelbergSpectralisOCT',eyeSide{eye});
        %find all vol files in directory
        allVols = dir(fullfile(inDirSub,'*.vol'));
        count = 0;
        
        %load all the data and display
        imIn = cell(1,3);
        figure(2)
        clf
        counter = 0;
        disp(['Subject #' num2str(s) ': ' subID ' ' eyeSide{eye}])
        for i = 1:length(allVols)
            PATH = fullfile(inDirSub,allVols(i).name);
            [HEADER, BSCANHEADER, SLO, Bscans] = openVolFast(PATH,'nodisp');
            
            %we're using this to filter out non-horizontal images remove if
            %you want to see all images
            if(HEADER.NumBScans > 1  || HEADER.ScanPattern > 1 || abs(BSCANHEADER.StartX - BSCANHEADER.EndX) < 1)
                continue
            else
                counter = counter+1;
            end
            Bscans = Bscans.^.25;
            subplot(2,ceil(length(allVols)/3),counter)
            imshow(Bscans(:,:,1));
            title(allVols(i).name, 'Interpreter', 'none')
            disp(allVols(i).name)
        end
        suptitle(['Subject #' num2str(s) ': ' subID ' ' eyeSide{eye}])

        pause;
        
    end
end
