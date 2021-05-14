clear all;close all;clc


load('D:\ycruzcav\tours\PET_unmixing\mat\initializations\brain_init_real_seg_02_02_3classes')
z1 = PET.Y;
load('D:\ycruzcav\ycruzcav_git\pet_unmixing\code\PET_unmixing_psf\mat\initializations\brain_init_real_5class')
z2 = PET.Y;

z1 = reshape(z1(15,:)',PET.W,PET.H,PET.D);
z1 = permute(z1,[2 1 3]);
z2 = reshape(z2(15,:)',PET.W,PET.H,PET.D);

figure;
% subplot(1,2,1)
% imagesc(z1(:,:,32))
% subplot(1,2,2)
% imagesc(z2(:,:,32))
plot3d(z1)
plot3d(z2)