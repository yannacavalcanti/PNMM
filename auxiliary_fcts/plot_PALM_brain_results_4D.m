function plot_PALM_brain_results_4D(PET,A_palm,B_palm,M_palm,A_gt,B_gt,M_gt,A0,M0,B0,B_depict,obj_fct,K,type)

%--------------------------------------------------------------
% Adjustements for plot
%--------------------------------------------------------------
% clear all
% load results\results_4D_psf_noise30_novar
A_palm = permute(reshape(A_palm',PET.W,PET.H,PET.D,K),[2 1 3 4]); %(PET.H|PET.W|K) : abundance cube
B_palm = permute(reshape(B_palm',PET.W,PET.H,PET.D,K-1),[2 1 3 4]); %(PET.H|PET.W|K) : abundance cube

A0 = permute(reshape(A0',PET.W,PET.H,PET.D,K),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
B0 = permute(reshape(B0',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
B_depict = permute(reshape(B_depict',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube

A_gt = permute(reshape(A_gt',PET.W,PET.H,PET.D,K),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
B_gt = permute(reshape(B_gt',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube

x_M = PET.time;

%% Plot

figure;
hold on
for i = 1:K-1
    subplot(K-1,4,1+4*(i-1));
    imagesc((B_palm(:,:,32,i)));
    title('Binding potential CNU')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    %             colormap gray
    colorbar
    subplot(K-1,4,2+4*(i-1));
    imagesc((B0(:,:,32,i)));
    title('Initial')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    colormap gray
    subplot(K-1,4,3+4*(i-1));
    imagesc((B_depict(:,:,32,i)));
    title('Initial')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    colormap gray
    subplot(K-1,4,4+4*(i-1));
    imagesc((B_gt(:,:,32,i)));
    title('Ground truth')
    caxis([0 max(abs(B_gt(:)))])
    %             caxis([-80 80])
    %             colormap gray
    colorbar
end
hold off

figure;
hold on
for i = 1:K
    subplot(K,4,1+(i-1)*4);
    imagesc(A_palm(:,:,32,i));
    title('PALM abundances')
    caxis([0 1]);
    %     colorbar
    %     colormap gray
    subplot(K,4,2+(i-1)*4);
    imagesc(A0(:,:,32,i));
    title('Initial')
    caxis([0 1]);
    %     colorbar
    %     colormap gray
    subplot(K,4,3+(i-1)*4);
    imagesc(A_gt(:,:,32,i));
    title('Ground truth')
    caxis([0 1]);
    %     colorbar
    %     colormap gray
    subplot(K,4,4+(i-1)*4);
    hold on
    plot(x_M,squeeze(M0(:,i)),'*k');
    plot(x_M,squeeze(M_palm(:,i)),'o-b');
    plot(x_M,squeeze(M_gt(:,i)),'+r');
    legend('Initial','PALM','Original')
    hold off
    title('Endmembers')
    axis([x_M(1) x_M(end) 0 max(M_palm(:))]);
    
end
hold off


figure;
subplot(1,3,1)
imagesc(squeeze(A_palm(64,:,:,3))')
title('PALM abundances - sagital')
subplot(1,3,2)
imagesc(squeeze(A0(64,:,:,3))')
title('Initial - sagital')
subplot(1,3,3)
imagesc(squeeze(A_gt(64,:,:,3))')
title('Ground truth - sagital')


figure;
plot(log(2:length(obj_fct(5,:))),log(obj_fct(5,2:end)))
title('Objective function')