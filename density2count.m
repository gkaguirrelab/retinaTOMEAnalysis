function output = density2count(image,R,numSectors,numRings)


% Terpin & McKendrick use:
% 15 segments per 90 deg of polar angle
% 2 rings per degree of ecc

for sec = 1:numSectors
    
    
    
    for ecc = 1:numRings
    binsize =     
    bin_pad= binsize/2;
    cx=size(I,1)/2;cy=size(I,2)/2;ix=size(I,1);iy=size(I,2);r=rad(i)+bin_pad;
    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    c_outer=((x.^2+y.^2)<=r^2);

    cx=size(I,1)/2;cy=size(I,2)/2;ix=size(I,1);iy=size(I,2);r=rad(i)-bin_pad;
    [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
    c_inner=((x.^2+y.^2)<=r^2);

    c_ring = c_outer-c_inner;   

        
        
    end
    
    
    
end
    