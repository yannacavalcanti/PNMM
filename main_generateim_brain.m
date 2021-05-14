clc;
close all;
clear all;

addpath generate_mixed_image_tac
addpath D:\ycruzcav\Database\mat
addpath plot\aux_fct plot\altmany-export_fig-5be2ca4
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load mat\generation_initialization\brain_initialization_4d

% Database size
DB = 1000;
% Number of eigenvectors
V = 3;
% SNR in dB
SNR =0;
% Dimensions
PET.W = 128;
PET.H = 128;
PET.D = 64;
% PSF: yes - 1, no - 0
PSF = 1;
PSF_type='3D FWHM';
% Noise type
noise_type='gaussian';% 'poisson','gaussian' or 'gaussian_cov'
Tracer_type = 'DPA'; % 'PE2I'
% Number of regions/classes
K=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------ Create synthetic image -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = 0;
% while(a PET.Y== 0)
% try
[PET,A_gt,M_gt,B_gt,B_basic,alpha_gt,E_gt,BP_gt] = generate_mixed_image_tac(PET,A0,M0,V,K,SNR,PSF,PSF_type,noise_type,Tracer_type);
%  a = 1;
% catch
%     a = 0;
% end
% end

%%% Find all TACs
% TACs = repmat(M_gt(:,1:end-1),1,1,V,2);
% for i=1:V
%     for k=1:K-1
%     TACs(:,k,i,:) = TACs(:,k,i,:)+(E_gt(:,:,k,i)*M_gt(:,1:end-1))*B_basic(k,:,i);
%     end
% end
% TACs = reshape(TACs,PET.L,(K-1)*V);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% t = unique(PET.Y','rows');
% figure(1)
% hold on;
% plot(PET.time,t','r')
% plot(PET.time,squeeze(M_gt),'b');
% hold off;
%%
lx=10; % width of subplot
ly=8; % height of subplot
margex=[1 1 1]*3; % margins around subplots
margey=[1 1 1]*3;
marge=0; % don't remember
orientation='landscape';
FS=24;
nc=1;
nr=1;
name='plot\images\TACs_exp';
[hh1, AX] = figaxes(nc, nr, lx, ly, margex, margey, marge, orientation,[0 0 1 1]);
%
figure(hh1);
axes(AX{1}{1});
hold on;
for k=1:2
    plot(PET.time,M_gt(:,k),'b','LineWidth',2)
    for n = 1:2
        t = M_gt(:,k)*(1+B_basic(k,n,1));
        for i=2:V
            t = t+(E_gt(:,:,k,i-1)*M_gt(:,k))*B_basic(k,n,i);
        end
        plot(PET.time,t,'r','LineWidth',2)
    end
end
plot(PET.time,M_gt(:,3),'b','LineWidth',2)
ylabel('Concentration(kBq/cc)','Interpreter','latex')
xlabel('Time(s)','Interpreter','latex')
set(gca,'XGrid','off','XMinorGrid','off','FontSize',FS,'Box','on','Layer','top');
   
axis([PET.time(1) PET.time(end) 0 max(M_gt(:))]);
title('Elementary TACs','Fontsize',FS,'FontName', 'Helvetica','Interpreter','latex');
hold off;



print('-depsc', '-r256',[name,'.eps']);
eps2pdf([name,'.eps'], [name,'.pdf'], 1, 0, 0, 100, {'-dSubsetFonts=true','-dEmbedAllFonts=true','-dCompatibilityLevel=1.4','-dPDFSETTUBGS=/printer'})
delete([name,'.eps'])


%%
lx=10; % width of subplot
ly=8; % height of subplot
margex=[1 1 1]*3; % margins around subplots
margey=[1 1 1]*3;
marge=0; % don't remember
orientation='landscape';
FS=24;
nc=1;
nr=1;
name='plot\images\TACs_exp_gray';
[hh1, AX] = figaxes(nc, nr, lx, ly, margex, margey, marge, orientation,[0 0 1 1]);
%
figure(hh1);
axes(AX{1}{1});
hold on;
for k=1
    plot(PET.time,M_gt(:,k),'b','LineWidth',2)
    for n = 1:2
        t = M_gt(:,k)*(1+B_basic(k,n,1));
        for i=2:V
            t = t+(E_gt(:,:,k,i-1)*M_gt(:,k))*B_basic(k,n,i);
        end
        plot(PET.time,t,'r','LineWidth',2)
    end
end
legend('nSB','SB')
ylabel('Concentration(kBq/cc)','Interpreter','latex')
xlabel('Time(s)','Interpreter','latex')
set(gca,'XGrid','off','XMinorGrid','off','FontSize',FS,'Box','on','Layer','top');
   
axis([PET.time(1) PET.time(end) 0 max(M_gt(:))]);
% title('Elementary TACs','Fontsize',FS,'FontName', 'Helvetica','Interpreter','latex');
hold off;



print('-depsc', '-r256',[name,'.eps']);
eps2pdf([name,'.eps'], [name,'.pdf'], 1, 0, 0, 100, {'-dSubsetFonts=true','-dEmbedAllFonts=true','-dCompatibilityLevel=1.4','-dPDFSETTUBGS=/printer'})
delete([name,'.eps'])

%%
lx=10; % width of subplot
ly=8; % height of subplot
margex=[1 1 1]*3; % margins around subplots
margey=[1 1 1]*3;
marge=0; % don't remember
orientation='landscape';
FS=24;
nc=1;
nr=1;
name='plot\images\TACs_exp_white';
[hh1, AX] = figaxes(nc, nr, lx, ly, margex, margey, marge, orientation,[0 0 1 1]);
%
figure(hh1);
axes(AX{1}{1});
hold on;
for k=2
    plot(PET.time,M_gt(:,k),'b','LineWidth',2)
    for n = 1:2
        t = M_gt(:,k)*(1+B_basic(k,n,1));
        for i=2:V
            t = t+(E_gt(:,:,k,i-1)*M_gt(:,k))*B_basic(k,n,i);
        end
        plot(PET.time,t,'r','LineWidth',2)
    end
end
legend('nSB','SB')
ylabel('Concentration(kBq/cc)','Interpreter','latex')
xlabel('Time(s)','Interpreter','latex')
set(gca,'XGrid','off','XMinorGrid','off','FontSize',FS,'Box','on','Layer','top');
   
axis([PET.time(1) PET.time(end) 0 max(M_gt(:))]);
% title('Elementary TACs','Fontsize',FS,'FontName', 'Helvetica','Interpreter','latex');
hold off;



print('-depsc', '-r256',[name,'.eps']);
eps2pdf([name,'.eps'], [name,'.pdf'], 1, 0, 0, 100, {'-dSubsetFonts=true','-dEmbedAllFonts=true','-dCompatibilityLevel=1.4','-dPDFSETTUBGS=/printer'})
delete([name,'.eps'])

%%
Y0 = reshape(PET.Y(15,:),PET.H,PET.W,PET.D);
figure(4)
subplot(1,3,1);
imagesc(Y0(:,:,32));
subplot(1,3,2);
imagesc(squeeze(Y0(:,64,:))');
subplot(1,3,3)
imagesc(squeeze(Y0(64,:,:))');

A0 = permute(reshape(A_gt',PET.W,PET.H,PET.D,K),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
figure;
hold on
for i = 1:K
    subplot(K,2,1+(i-1)*2);
    imagesc(A0(:,:,32,i));
    title('Abundances')
    caxis([0 1])
    colorbar
    subplot(K,2,2+(i-1)*2);
    hold on
    plot(PET.time,squeeze(M_gt(:,i)),'o-b');
    %     legend('Initial')
    hold off
    title('Endmembers')
    %     axis([0 20 min(M_gt(:)) max(M_gt(:))])
    
end
hold off

BP0 = permute(reshape(BP_gt',PET.W,PET.H,PET.D,K-1),[2 1 3 4]);           %(PET.H|PET.W|K) : abundance cube
figure;
hold on
for i = 1:K-1
    subplot(K-1,1,i);
    imagesc(BP0(:,:,32,i));
    title('Binding potential')
    %     caxis([-max(abs(B(:))) max(abs(B(:)))])
    %     colorbar hot
end
hold off
%% B
lx=10; % width of subplot
ly=8; % height of subplot
margex=[1 1 1]*3; % margins around subplots
margey=[1 1 1]*3;
marge=0; % don't remember
orientation='landscape';
FS=24;
nc=K-1;
nr=1;
name='images\BP_ground_truth';
[hh1, AX_BP] = figaxes(nc, nr, lx, ly, margex, margey, marge, orientation,[0 0 1 1]);
%
figure(hh1);
tissueBP = {'Gray matter','White matter'};
for i=1:K-1
axes(AX_BP{i}{1});
h = imagesc(BP0(12:71,21:100,30,i)');
set(h,'alphadata',~isnan(BP0(12:71,21:100,30,i)'))
colorbar;
set(gca, 'ydir', 'reverse');
axis tight
title(tissueBP{i},'Fontsize',FS,'FontName', 'Helvetica','Interpreter','latex');
caxis([0 max(BP0(:))]);
axis off
colormap default
set(gca,'XGrid','off','XMinorGrid','off','FontSize',FS,'Box','on','Layer','top');
end
print('-depsc', '-r256',[name,'.eps']);
eps2pdf([name,'.eps'], [name,'.pdf'], 1, 0, 0, 100, {'-dSubsetFonts=true','-dEmbedAllFonts=true','-dCompatibilityLevel=1.4','-dPDFSETTUBGS=/printer'})
delete([name,'.eps'])



PET.C = max(PET.Y(:));

clear A0 B0 DB deltaM i K PET.L m_dynamic Nv m_var SNR_ROI a M0 A M AB
save(strcat('mat\generated_images\brain_synt_4d_exp_',Tracer_type,'_',noise_type,'_noise',int2str(SNR)),'PET','A_gt','M_gt','B_gt','B_basic','alpha_gt','E_gt','BP_gt')
