function moving_reg = applyTransformation(moving,target,T_bck)

ur = 1:size(target,2);
vr = 1:size(target,1);
[u,v] = meshgrid(ur,vr) ;
H = T_bck;
z_ = H(3,1) * u + H(3,2) * v + H(3,3);
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
moving_reg = vl_imwbackward(im2double(moving),u_,v_);