function nii = convertVoltoNii(inputVolPath,outputDirPath)



[header, BScanHeader, slo, BScans] = openVolFast (inputVolPath);

origin = [0 0 0];
img = BScans.^(.25);%we cube root intensities to allow us to see the contrast
img(img>5) = 0;
img = rot90(img,-1);%Rotate to reorient the image. Remember, this flips the index for the depth direction.
res = [header.ScaleX header.ScaleZ header.Distance];
nii = make_nii(img, res, origin);

[filepath,name,ext] = fileparts(inputVolPath);
save_nii(nii, fullfile(outputDirPath,[name '.nii']))
