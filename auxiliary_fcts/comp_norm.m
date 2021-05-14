function norm_HHt = comp_norm(H,W)
% Compute smoothing norm with spatial penalization
% 
%%%%%% Inputs:
% - W,H    width and height of image
%%%%%% Outputs:
% - norm_HHt    the norm of HHt
%
% Yanna Cavalcanti, April 2016

HHt_hor = 2*((W-2)*36+2*(W-1)*4+2*4)+(H-2)*((W-2)*64+2*(W-1)*4+2*36);
HHt_ver = 2*(H-1)*(W*4);


norm_HHt = sqrt(HHt_ver+HHt_hor);