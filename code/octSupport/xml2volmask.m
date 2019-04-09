function mask = xml2volmask(xmlfile)
% convert xml from OCT Explorer to volume mask
% syntax: mask = xml2volmask(xmlfile)
%	input: xmlfile - xml filepath
% 	output: mask - volume mask of each slab (1 to 10) (3D matrix: nx*ny*nz)
% written by Jin Gahm, LONI, USC

fin = fopen(xmlfile,'r');
l = 0;
while ~feof(fin)
    
    ss = fgetl(fin);
    
    if(strfind(ss,'<size>'))
        fgetl(fin);
        s = fgetl(fin);
        ny = str2double(s((strfind(s,'x')+2):(strfind(s,'/')-2)));
        s = fgetl(fin);
        nz = str2double(s((strfind(s,'y')+2):(strfind(s,'/')-2)));
        s = fgetl(fin);
        nx = str2double(s((strfind(s,'z')+2):(strfind(s,'/')-2)));
    end
    
    if(strfind(ss,'<surface_num>'))
        r = ss((strfind(ss,'>')+1):(strfind(ss,'/')-2));
        nlayers = str2double(r);
        bd_pts = zeros(nx,ny,nlayers);
        undefined_regions = false(nx,ny);
    end
    
    if(strfind(ss,'label'))
        fgetl(fin);
        fgetl(fin);
        l = l+1;
        for i=0:nx-1
            fgetl(fin);
            for j=0:ny-1
                s = fgetl(fin);
                r = s((strfind(s,'y')+2):(strfind(s,'/')-2));
                %fprintf(fout,'%d %d %s\n',i,j,r);
                bd_pts(i+1,j+1,l) = str2double(r)+1;
            end
            fgetl(fin);
        end
    end
    
    if(strfind(ss,'<undefined_region>'))
        while (strfind(fgetl(fin),'<ascan>'))
            s = fgetl(fin);
            x = str2double (s((strfind(s,'x')+2):(strfind(s,'/')-2)))+1;
            s = fgetl(fin);
            z = str2double(s((strfind(s,'z')+2):(strfind(s,'/')-2)))+1;
            undefined_regions(z,x) = 1;
            fgetl(fin);
        end
    end
end
fclose(fin);

%%
bd_pts = permute(bd_pts,[2 1 3]);
nx = size(bd_pts,1);
ny = size(bd_pts,2);
nlayers = size(bd_pts,3);
mask = uint8(zeros(nx,ny,nz));

for i=1:nx
    for l=1:nlayers-1
        x1 = false(ny,nz);
        x2 = false(ny,nz);
        p1 = bd_pts(i,:,l);
        p2 = bd_pts(i,:,l+1);
        for j=1:ny-1
            q1 = line2d([j round(p1(j))],[j+1 round(p1(j+1))]);
            x1(sub2ind([ny nz],q1(:,1),q1(:,2))) = 1;
            q2 = line2d([j round(p2(j))],[j+1 round(p2(j+1))]);
            x2(sub2ind([ny nz],q2(:,1),q2(:,2))) = 1;
        end
        
        x3 = true(ny+2,nz);
        x3(2:end-1,:) = x1+x2;
        x4 = imfill(x3,'holes');
        x5 = x4(2:end-1,:);
        j = find(x5);
        mask(i,j) = l*x5(j);
    end
end

mask = mask.*repmat(uint8(~undefined_regions'),[1 1 nz]);

%mask = permute(mask,[2 1 3]);
mask = flipdim(mask,3);
mask = flipdim(mask,1);
mask = flipdim(mask,2);

end % xml2volmask

%% Local function
function p = line2d(p1,p2)

nPixels=max(abs(p1-p2))+1;
X=linspace(p1(2), p2(2), nPixels);
f=(p1(1)-p2(1))/(p1(2)-p2(2));
if isinf(f)
    Y=linspace(p1(1),p2(1),nPixels);
else
    f(2)=-det([p1;p2])/(p1(2) - p2(2));
    Y=(f(1)*X+f(2));
end
p = [round(Y);round(X)]';
end