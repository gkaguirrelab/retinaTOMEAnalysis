function [T_bck, moving_reg, numMatches_all]  = siftReg(moving,target,TransType,saveFlag,outputDir)
%Written by Min Chen (minchen1@upenn.edu)

MN = 1;%using 3 modalities in this example
N = 2;%two images to match


%save SIFT features
f_all = cell(MN,N);
d_all = cell(MN,N);

%read images
im1 = cell(MN,1);
im2 = cell(MN,1);
im1{1}=target;
im2{1}=moving;

%find SIFT features
for m=1:MN
    im = im2single(im1{m});
    [f1,d1] = vl_sift(im,'Levels',55);
    [f1, d1] = filterSiftFeaturesByROI(im, f1, d1, 0);
    f_all{m,1} = f1;
    d_all{m,1} = d1;
    
    im = im2single(im2{m});
    [f2, d2] = vl_sift(im,'Levels',55);
    [f2, d2] = filterSiftFeaturesByROI(im, f2, d2, 0);
    f_all{m,2} = f2;
    d_all{m,2} = d2;
end
%saveFlag = 0;
if(saveFlag)
    mkdir(outputDir);
end
[T_bck, numOkMatches_all, numMatches_all]=sift_mosaic_fast_MultiModal(im1, im2, outputDir,saveFlag,f_all(:,1),d_all(:,1),f_all(:,2),d_all(:,2),TransType);
%apply transformation 
moving_reg = applyTransformation(moving,target,T_bck);

