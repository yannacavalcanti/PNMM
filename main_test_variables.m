clc,
clear all, format compact;
close all
addpath palms auxiliary_fcts;
addpath D:\ycruzcav\Database\mat
addpath nonlinear_PETunmixing
addpath mat

% --------------------------------------------------------------
% Choose initialization
% --------------------------------------------------------------
tracer_type = 'DPA';
img_name = strcat('_4d_exp_',tracer_type,'_gaussian_noise0');
init_type = 'segmentation'; % type of initialization: 'nfindr','vca','ROI nfindr','ROI vca','rand','zero','ground truth','test A','test B','test M'
type_init = 'create'; % 'create','load'
img_type='synt';
K = 3;
V = 3;

% --------------------------------------------------------------
% Choose image
% --------------------------------------------------------------

[PET,A_gt,B_gt,M_gt,alpha_gt,BP_gt,B_basic,vect_prop] = choose_image(img_type,img_name,K);

%--------------------------------------------------------------
% Unmixing parameters
%--------------------------------------------------------------

% PALM parameters
Niter = 3000;
beta = 0.1;   % endmember penalization parameter: 0.0086 - pen: 4.37, fct_term = 0.0916
lbd = 0.5; % variability penalization - sparsity
gamma= 0;
eta = 0.5; %abundance penalization - smoothing
alpha_lim = [0.6,0.05];
B_lim = [eps,-0.25,eps;0.7,-eps,0.15];
epsilon =1e-03;
C = max(var(PET.Y(:,PET.mask)));
%--------------------------------------------------------------
% Adjustments
%--------------------------------------------------------------

A0 = A_gt;
M0 = M_gt;
B0 = B_gt;
alpha0 = alpha_gt;

% A0(:,PET.mask) = rand(size(A_gt(:,PET.mask)));
M0 = zeros(size(M_gt));
% B0(:,PET.mask,1) = eps*ones(size(B_gt(:,PET.mask,1)));
% B0(:,PET.mask,2) = -eps*ones(size(B_gt(:,PET.mask,2)));
% B0(:,PET.mask,3) = eps*ones(size(B_gt(:,PET.mask,3)));
% alpha0 = [0.6*ones(1,2);0.01*ones(1,2)];

BP0 = B0(:,:,1);
for k=1:K-1
    BP0(k,:) = BP0(k,:)+sum(squeeze(B0(k,:,2:end))./repmat(alpha0(:,k)',size(A0,2),1),2)';
end
%--------------------------------------------------------------
% PNMM
%--------------------------------------------------------------
[A_cnu,M_cnu,B_cnu,alpha_cnu,BP_cnu,Q,elapsedTime,obj_fct,ni] = nonlin_PETunmixing(PET,A0,M0,B0,alpha0,Niter,lbd*C,gamma*C,beta*C,eta*C,alpha_lim,B_lim,epsilon);
%--------------------------------------------------------------
% Compute SNR
%--------------------------------------------------------------
N_it = ni-1;
load('results\results_4D_GT')
% [err0,err0_MA] = cmpt_all_err_norm(PET,A0,M0,B0(:,:,1),alpha0,BP0,A_gt,M_gt,B_gt(:,:,1),alpha_gt,BP_gt,K,V);
% [err_cnu,err_cnu_MA] = cmpt_all_err_norm(PET,A_cnu,M_cnu,B_cnu(:,:,1),alpha_cnu,BP_cnu,A_gt,M_gt,B_gt(:,:,1),alpha_gt,BP_gt,K,V);

% --------------------------------------------------------------
% Plot and save
%--------------------------------------------------------------


plot_PALM_brain_results_4D(PET,A_cnu,BP_cnu,M_cnu,A_gt,BP_gt,M_gt,A0,M0,BP0,BP0,obj_fct,K,img_type);
save(strcat('results\results_4D_',img_type,img_name))
% 1e-2
% err_cnu =
%     0.0489    0.0815    0.4226    0.0536
% 5e-3
% err_cnu =
%     0.0897    0.0858    0.4037    0.0803
