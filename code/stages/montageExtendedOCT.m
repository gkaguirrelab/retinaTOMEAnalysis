function montageExtendedOCT(dataDir,scanInfoFile,saveDir)
%Purpose: Montages the extended horizontal OCT into a single aligned image
%Inputs:
%dataDir - Input directory with all OCT data
%scanInfoFile - spreadsheet describing which scans are the
%left/center/right images in the montage
%saveDir - output directory to save the montaged images

[NUM,scanInfoTXT,RAW]=xlsread(scanInfoFile);

FullSave = 0;%flat for saving additional information about the alignments
allSubs = dir(fullfile(dataDir,'1*'));%find all subject directories
eyeSide = {'OS', 'OD'};
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
        
        
        %first we load the data using the spreadsheet information
        imIn = cell(1,3);

        if (eye == 1)%OS
            subVols = scanInfoTXT(s+2,3:5);
        else%OD
            subVols = scanInfoTXT(s+2,7:9);
        end
        
        %load all three images
        count = 0;
        for i = 1:3
            count = count+1;
            PATH = fullfile(inDirSub,subVols{i});
            [HEADER, BSCANHEADER, SLO, Bscans] = openVolFast(PATH);
            Bscans = Bscans.^.25;
            subplot(3,1,count);            
            imshow(Bscans(:,:,1));
            imIn{count} = Bscans(:,:,1);
            res = [HEADER.ScaleZ HEADER.ScaleX .001];
        end
        
      
        %set images
        im1=imIn{1};
        im2=imIn{2};
        im3=imIn{3};
        
        %register to center image
        [T_bck_im2_im1, im2_reg, numMatches_all]  = siftReg(im2,im1,1,FullSave,fullfile(saveDirSub,'2-1'));
        [T_bck_im3_im1, im3_reg, numMatches_all]  = siftReg(im3,im1,1,FullSave,fullfile(saveDirSub,'3-1'));
        
        %combine all 3 images together        
        %first find the min/max bounds across all 3 images
        H=T_bck_im2_im1;
        box2 = [1  size(im2,2) size(im2,2)  1              ;
            1  1              size(im2,1)  size(im2,1) ;
            1  1              1               1            ] ;
        box2_ = inv(H) * box2 ;
        box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
        box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
        
        H=T_bck_im3_im1;
        box3 = [1  size(im3,2) size(im3,2)  1              ;
            1  1              size(im3,1)  size(im3,1) ;
            1  1              1               1            ] ;
        box3_ = inv(H) * box3 ;
        box3_(1,:) = box3_(1,:) ./ box3_(3,:) ;
        box3_(2,:) = box3_(2,:) ./ box3_(3,:) ;
        
        %create the final image grid for the montage
        ur = min([1 box2_(1,:) box3_(1,:)]):max([size(im1,2) box2_(1,:) box3_(1,:)]) ;
        vr = min([1 box2_(2,:) box3_(2,:)]):max([size(im1,1) box2_(2,:) box3_(2,:)]) ;
        [u,v] = meshgrid(ur,vr) ;
        
        %place im1 in combined space
        im1_ = vl_imwbackward(im2double(im1),u,v) ;
        
        %place registered im2 in combined space
        H = T_bck_im2_im1;
        z_ = H(3,1) * u + H(3,2) * v + H(3,3);
        u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
        v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
        im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
        im2_ = im2_(:,:,1);
        
        %place registered im3 in combined space
        H = T_bck_im3_im1;
        z_ = H(3,1) * u + H(3,2) * v + H(3,3);
        u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
        v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
        im3_ = vl_imwbackward(im2double(im3),u_,v_) ;
        im3_ = im3_(:,:,1);
        
        %average all 3 images
        overlap= ((im1_ ~= 0) & (~isnan(im1_))) + ((im2_ ~= 0) & (~isnan(im2_))) + ((im3_ ~= 0) & (~isnan(im3_)));
        overlap(overlap == 0) = 1;
        
        %remove nan's in the image
        im1_(isnan(im1_))=0;
        im2_(isnan(im2_))=0;
        im3_(isnan(im3_))=0;
        imAvg = (im1_+im2_+im3_)./(overlap);
        imOut = cat(3,im1_, im2_,im3_);
        
        %save out images
        if(FullSave)%this saves out additional registation information
            saveTif(im1,saveDirSub,'im1.tif');
            saveTif(im2,saveDirSub,'im2.tif');
            saveTif(im3,saveDirSub,'im3.tif');
            
            saveTif(im1_,saveDirSub,'im1_Trans.tif');
            saveTif(im2_,saveDirSub,'im2_Trans.tif');
            saveTif(im3_,saveDirSub,'im3_Trans.tif');
        end
        %Save montaged average
        saveTif(imAvg,saveDirSub,'AvgAll_Trans.tif');

        %save .mat with transformation information
        save(fullfile(saveDirSub,[subID '_' eyeSide{eye} '_ExtendedOCT.mat']),'imIn','T_bck_im2_im1','T_bck_im3_im1','imAvg','imOut','subID','eye','res');

        %export as nii for delineation
        origin = [0 0 0];
        
        %Reorient matlab img to nii img
        img = permute(imOut,[2 1 3]);
        img = rot90(img,2);
        niires = [res(2) res(1) res(3)];
        nii = make_nii(img, niires, origin);
        save_nii(nii, fullfile(saveDirSub,[subID '_' eyeSide{eye} '_ExtendedOCT.nii']))
    end
end