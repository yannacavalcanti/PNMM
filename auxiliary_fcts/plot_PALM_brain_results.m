function plot_PALM_brain_results(PET,A_palm,B_palm,M_palm,A_gt,B_gt,M_gt,A0,M0,B0,obj_fct,err_A,err_M,err_beta,K,Nv,type)
if strcmp(type,'synt')
%--------------------------------------------------------------
% Adjustements for plot
%--------------------------------------------------------------

A_palm = permute(reshape(A_palm',PET.W,PET.H,K),[2 1 3]); %(PET.H|PET.W|K) : abundance cube    
B_palm = permute(reshape(B_palm',PET.W,PET.H,Nv),[2 1 3]); %(PET.H|PET.W|K) : abundance cube 
          
A0 = permute(reshape(A0',PET.W,PET.H,K),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    
B0 = permute(reshape(B0',PET.W,PET.H,Nv),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    

A_gt = permute(reshape(A_gt',PET.W,PET.H,K),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    
B_gt = permute(reshape(B_gt',PET.W,PET.H,Nv),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    

vect_frame_time=[60 60 60 60 60 120 120 120 120 120 150 150 300 300 300 300 300 300 300 300];
x_M = cumsum(vect_frame_time);

%% Plot

figure;
hold on
for i = 1:Nv
    subplot(Nv,3,1+(i-1)*2);
    imagesc((B_palm(:,:,i)));
    title('PALM internal abundances')
%     caxis([0 max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    caxis([-max(max(abs(B_palm(:))),max(abs(B_gt(:)))) max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    colormap hot
    colorbar
    subplot(Nv,3,2+(i-1)*2);
    imagesc((B0(:,:,i)));
    title('Initial')
    caxis([-max(max(abs(B_palm(:))),max(abs(B_gt(:)))) max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    colormap gray
    subplot(Nv,3,3+(i-1)*2);
    imagesc((B_gt(:,:,i)));
    title('Ground truth')
%     caxis([0 max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    caxis([-max(max(abs(B_palm(:))),max(abs(B_gt(:)))) max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    colormap hot
    colorbar
end
hold off

figure;
hold on
for i = 1:K
    subplot(K,4,1+(i-1)*4);
    imagesc(A_palm(:,:,i));
    title('PALM abundances')
    caxis([0 1])
    %     colorbar
%     colormap gray
    subplot(K,4,2+(i-1)*4);
    imagesc(A0(:,:,i));
    title('Initial')
    caxis([0 1])
    %     colorbar
    subplot(K,4,3+(i-1)*4);
    imagesc(A_gt(:,:,i));
    title('Grount truth')
    caxis([0 1])
    %     colorbar
%     colormap gray
    subplot(K,4,4+(i-1)*4);
    hold on
    plot(x_M,squeeze(M0(:,i)),'*k');
    plot(x_M,squeeze(M_palm(:,i)),'o-b');
    plot(x_M,squeeze(M_gt(:,i)),'+r');
    legend('Initial','PALM','Ground truth')
    hold off
    title('Endmembers')
    axis([x_M(1) x_M(end) 0 1.3])

end
hold off

figure;
subplot(2,2,1)
plot(log10(2:size(obj_fct,2)),log10(obj_fct(5,2:end))')
legend('Data fitting term','A term','M term','B term','Objective fct')
title('Objective function')
subplot(2,2,2)
plot(log10(2:length(err_A)),log10(err_A(2:end)))
title('Abundance error')
subplot(2,2,3)
plot(log10(2:length(err_M)),log10(err_M(2:end)))
title('Endmember error')
subplot(2,2,4)
plot(log10(2:length(err_beta)),log10(err_beta(2:end)))
title('Variability error')
% 
% figure;
% subplot(1,3,1)
% plot(log(2:length(obj_fct)),log(obj_fct(2:end)))
% title('Objective function')
% subplot(1,3,2)
% plot(log(2:length(err_A)),log(err_A(2:end)))
% title('Abundance error')
% subplot(1,3,3)
% plot(log(2:length(err_M)),log(err_M(2:end)))
% title('Endmember error')
else
    %--------------------------------------------------------------
% Adjustements for plot
%--------------------------------------------------------------

A_palm = permute(reshape(A_palm',PET.W,PET.H,K),[2 1 3]); %(PET.H|PET.W|K) : abundance cube    
B_palm = permute(reshape(B_palm',PET.W,PET.H,Nv),[2 1 3]); %(PET.H|PET.W|K) : abundance cube 
          
A0 = permute(reshape(A0',PET.W,PET.H,K),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    
B0 = permute(reshape(B0',PET.W,PET.H,Nv),[2 1 3]);           %(PET.H|PET.W|K) : abundance cube    

vect_frame_time=[60 60 60 60 60 120 120 120 120 120 150 150 300 300 300 300 300 300 300 300];
x_M = cumsum(vect_frame_time);

%% Plot

figure;
hold on
for i = 1:Nv
    subplot(Nv,2,1+(i-1)*2);
    imagesc(abs(B_palm(:,:,i)));
    title('PALM internal abundances')
    caxis([0 max(max(abs(B_palm(:))),max(abs(B_gt(:))))])
    colormap gray
    colorbar
    subplot(Nv,2,2+(i-1)*2);
    imagesc(abs(B0(:,:,i)));
    title('Initial')
    caxis([0 max(max(abs(B_palm(:))),max(abs(B0(:))))])
    colormap gray
end
hold off

figure;
hold on
for i = 1:K
    subplot(K,3,1+(i-1)*3);
    imagesc(A_palm(:,:,i));
    title('PALM abundances')
    caxis([0 1])
    %     colorbar
%     colormap gray
    subplot(K,3,2+(i-1)*3);
    imagesc(A0(:,:,i));
    title('Initial')
    caxis([0 1])
    %     colorbar
    subplot(K,3,3+(i-1)*3);
    hold on
%     plot(x_M,squeeze(M0(:,i)),'*k');
%     plot(x_M,squeeze(M_palm(:,i)),'o-b');
    legend('Initial','PALM')
    hold off
    title('Endmembers')
    axis([x_M(1) x_M(end) 0 1.3])

end
hold off

figure;
plot(log10(2:size(obj_fct,2)),log10(obj_fct(5,2:end))')
title('Objective function')


figure;
hold on
subplot(1,4,1);
imagesc(A_palm(:,:,1));
title('PALM abundances')
caxis([0 1])
%     colorbar
%     colormap gray
subplot(1,4,2);
imagesc(A0(:,:,1));
title('Initial')
caxis([0 1])
%     colorbar
subplot(1,4,2);
imagesc(A_gt(:,:,1));
title('Initial')
caxis([0 1])
%     colorbar
subplot(1,4,3);
hold on
plot(squeeze(M0(:,1)),'*k');
plot(squeeze(M_palm(:,1)),'o-b');
legend('Initial','PALM')
hold off
title('Endmembers')
axis([x_M(1) x_M(end) 0 1.3])

hold off

end
