clc;
close all;
clear all;

addpath generate_mixed_image_tac
addpath D:\ycruzcav\Database\mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load mat\generation_initialization\brain_initialization_4d
% Number of f rames
PET.L=size(M0,1);
% Number of eigenvectors
V = 3;
% Dimensions
PET.W = 128;
PET.H = 128;
PET.D = 64;
% Noise type
noise_type='gaussian_mult';% 'poisson','gaussian' or 'gaussian_cov'
Tracer_type = 'DPA'; % 'PE2I'
% Number of regions/classes
K=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------ Create synthetic image -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load data/dpa_tacs
vect_time = time_dpa;
PET.L=length(vect_time);
[gt_phantom,A,M,B,B_basic,alpha,E,BP] = generate_parameters_dpa(A0,M0,V,vect_time,K,M_dpa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbr = find(sum(BP)>0.3);
z = gt_phantom(:,sbr);
figure(3)
hold on
plot(vect_time,z(:,1:50:end),'b');
plot(vect_time,squeeze(M(:,1)),'r');
title('Specific Binding TACs')
xlabel('Time')
ylabel('Concentration')
hold off


gt_phantom = reshape(gt_phantom',PET.W,PET.H,PET.D,PET.L);
figure(4)
subplot(1,3,1);
imagesc(gt_phantom(:,:,32,15));
subplot(1,3,2);
imagesc(squeeze(gt_phantom(:,64,:,15))');
subplot(1,3,3)
imagesc(squeeze(gt_phantom(64,:,:,15))');


% keyboard
%clear A0 B0 DB deltaM m_dynamic SNR_ROI a M0 A M AB z Y0 V_type PSF PSF_type caso
A = reshape(A',PET.W,PET.H,PET.D,K);
BP = reshape(BP',PET.W,PET.H,PET.D,K-1);

% Find labels and TACs
z = reshape(gt_phantom,128*128*64,PET.L)';
z = z(:,sum(z)>0);
[idx,Tacs] = kmeans(z',7);
label_matrix = zeros(size(gt_phantom(:,:,:,15)));
label_matrix(gt_phantom(:,:,:,15)>0)=idx;

save('mat\generated_images\gt_phantom_PNMM','gt_phantom','vect_time','label_matrix','Tacs','A','M','B','B_basic','alpha','BP')
