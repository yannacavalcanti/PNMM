function [PET,A0,M0,M_ROI,vect_prop] = load_real_image(K)

load('C:\Users\yccru\Documents\yanna_phd\tours\PET_unmixing\mat\initializations\brain_init_real_seg_02_03_3classes','A0','PET','M0')

abund = reshape(permute(reshape(A0(1:2,:)',PET.W,PET.H,PET.D,K-1),[2 1 3 4]),PET.W*PET.H*PET.D,K-1)';
endm = M0(:,1:2);

load(strcat('C:\Users\yccru\Documents\yanna_phd\ycruzcav_git\pet_unmixing\code\PET_unmixing_psf\mat\initializations\brain_init_real_5class'),'A0','A_gt','B0','B_gt','M0','M_gt','PET','M_ROI','vect_prop')
abund2 = A0;
endm2 = M0;
mask2 = find(A0(5,:)==1);

z = zeros(1,size(PET.Y,2));
z(PET.mask)=1;
z(mask2)=0;
PET.mask = find(z==1);
K=3;

A0 = nan(K,PET.W*PET.H*PET.D);
A0(:,PET.mask)=0;
A0(1,(abund2(1,:)==1 & abund(2,:)==0) | (abund(1,:)==1) & abund2(3,:)==0)=1;
A0(2,(abund2(2,:)==1 & abund(1,:)==0) | (abund(2,:)==1) & abund2(3,:)==0)=1;
% A0(1:2,:)=abund;
A0(3,:)=abund2(3,:);

M0(:,1)=endm2(:,2);
M0(:,2)=endm2(:,4);
M0(:,3)=endm2(:,3);
vect_frame_time=[10*ones(1,6) 30*ones(1,8) 60*ones(1,4) 120*ones(1,5) 300*ones(1,8)];
PET.time = cumsum(vect_frame_time);