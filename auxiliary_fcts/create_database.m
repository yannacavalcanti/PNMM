
load D:\ycruzcav\Database\mat\ground_truth
load mat\brain_mask

Y = reshape(img_gt_ds,128*128*64,20)';
K=4;
A0 = zeros(K,size(Y,2));
Y = Y(:,V_mask);

%         Endmember initialization
[M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm(Y,K,'nfindr');
% Abundance initialization (SUNSAL [2])
A_Mask = sunsal(M0,Y,'POSITIVITY','yes','ADDONE','yes');
N = size(A_Mask,2);
% % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
for n = 1:N
    A_Mask(:,n) = ProjectOntoSimplex(A_Mask(:,n),1);
end
A0(:,V_mask) = A_Mask;

W=128;
H=128;
D = 64;

A0 = reshape(A0',W,H,D,K);           %(H|W|K) : abundance cube   

figure;
hold on;
for i = 1:K
colormap gray
    subplot(K,2,1+(i-1)*2);
    imagesc(A0(:,:,32,i));
    title('Initial')
    caxis([0 1])
%     colorbar
colormap gray
    subplot(K,2,2+(i-1)*2);
    plot(M0(:,i))
end
hold off;

A0 = reshape(A0,W*H*D,K)';

abund = A0;
endm = M0;
abund(1,:) = A0(2,:);
abund(2,:) = A0(1,:);
endm(:,1) = M0(:,2);
endm(:,2) = M0(:,1);
A0 = abund;
M0 = endm;