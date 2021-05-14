clc;
clear all;
close all;
load 'D:\ycruzcav\Database\mat\pet_ground_truth'
addpath auxiliary_fcts mat
load brain_mask

K = 4;
Z = Y(:,V_mask);
A0 = zeros(K,128*128*64);
%         Endmember initialization
[M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm(Z,K,'nfindr');
plot(M0)
% Abundance initialization (SUNSAL [2])
A_mask= sunsal(M0,Z,'POSITIVITY','yes','ADDONE','yes');
N = size(A0,2);
% % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
for n = 1:length(V_mask)
    A_mask(:,n) = ProjectOntoSimplex(A_mask(:,n),1);
end
A0(:,V_mask) = A_mask;
keyboard
%
%         M = M0;
%         M(:,1) = M0(:,4);
%         M(:,4) = M0(:,1); plot(M)
%         M0 = M;


% 
% clear
% 

K=K+2;
     
A0 = reshape(A0,K,PET.W,PET.H,PET.D);           %(PET.H|PET.W|K) : abundance cube    
A_gt = reshape(A_gt,K,PET.W,PET.H,PET.D);           %(PET.H|PET.W|K) : abundance cube    

figure;
hold on
for i = 1:K
    subplot(K,2,1+(i-1)*2);
    imagesc(squeeze(A0(i,:,:,32)));
    title('Initial')
    %caxis([0 1]);
    subplot(K,2,2+(i-1)*2);
%     plot(M0(i,:))
    imagesc(squeeze(A_gt(i,:,:,32)));
    title('Ground truth')
    caxis([0 1]);
end
hold off
    
A0 = reshape(A0,K,PET.W*PET.H*PET.D);           %(PET.H|PET.W|K) : abundance cube    
A_gt = reshape(A_gt,K,PET.W*PET.H*PET.D);           %(PET.H|PET.W|K) : abundance cube    

abund = A0;
endm = M0;
abund(1,:) = A0(3,:);
% abund(2,:) = A0(1,:);
abund(3,:) = A0(4,:);
abund(4,:) = A0(1,:);
endm(:,1) = M0(:,3);
% endm(:,2) = M0(:,4);
endm(:,3) = M0(:,4);
endm(:,4) = M0(:,1);
A0 = abund;
M0 = endm;
% 
% abund = A_palm;
% endm = M_palm;
% abund(2,:) = A_palm(4,:);
% % % abund(3,:) = A0(4,:);
% abund(4,:) = A_palm(2,:);
% endm(:,2) = M_palm(:,4);
% % % endm(:,3) = M0(:,4);
% endm(:,4) = M_palm(:,2);
% A_palm = abund;
% M_palm = endm;
% 
% ind = 1:K;
% n = zeros(1,K);
% abund = zeros(size(A0));
% for i=1:K
%     for j = 1:length(ind)
%     n(j) = norm(M0(:,i)-M_gt(:,ind(j)));
%     end
%     im=find(min(n));
%     abund(i,:) = A0(ind(im),:);
%     ind(im) = [];
%     clear n
% end
% 


BP0 = permute(reshape(BP0',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
BP_gt = permute(reshape(BP_gt',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube

figure;
hold on
for i = 1:K-1
    subplot(K-1,3,2+3*(i-1));
    imagesc((BP0(:,:,32,i)));
    title('Initial')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    colormap gray
    subplot(K-1,3,3+3*(i-1));
    imagesc((BP_gt(:,:,32,i)));
    title('Ground truth')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    %             colormap gray
    colorbar
end
BP0 = reshape(BP0,K-1,PET.W*PET.H*PET.D);           %(PET.H|PET.W|K) : abundance cube    
BP_gt = reshape(BP_gt,K-1,PET.W*PET.H*PET.D);           %(PET.H|PET.W|K) : abundance cube    
