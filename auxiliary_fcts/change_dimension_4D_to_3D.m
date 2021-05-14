function [PET,A0,A_gt,B0,B_gt,B_mat] = change_dimension_4D_to_3D(PET,A0,A_gt,B0,B_gt,B_mat);

K = size(A0,1);
Nv = size(B0,1);
% Mask
I_Mask = zeros(1,size(PET.Y,2));
I_Mask(PET.mask) = 1;
I_Mask = reshape(I_Mask,PET.W,PET.H,PET.D);
I_Mask = I_Mask(:,:,32);
I_Mask = reshape(I_Mask,1,PET.W*PET.H);
PET.mask = find(I_Mask == 1);
% Mask
% I_Mask = zeros(1,size(PET.Y,2));
% I_Mask(PET.mask_nan) = 1;
% I_Mask = reshape(I_Mask,PET.W,PET.H,PET.D);
% I_Mask = I_Mask(:,:,32);
% I_Mask = reshape(I_Mask,1,PET.W*PET.H);
% PET.mask_nan = find(I_Mask == 1);

% Abundance
[A0] = four_to_threeD(A0,K,PET.W,PET.H,PET.D);
[A_gt] = four_to_threeD(A_gt,K,PET.W,PET.H,PET.D);

% % A*B
% [AB_gt] = four_to_threeD(AB_gt,Nv,PET.W,PET.H,PET.D);

% Variability
[B0] = four_to_threeD(B0,Nv,PET.W,PET.H,PET.D);
[B_gt] = four_to_threeD(B_gt,Nv,PET.W,PET.H,PET.D);
[B_mat] = four_to_threeD(B_mat,PET.L,PET.W,PET.H,PET.D);

% Image
[PET.Y] = four_to_threeD(PET.Y,PET.L,PET.W,PET.H,PET.D);
PET.D = 1;