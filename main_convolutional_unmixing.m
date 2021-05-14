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
img_name = strcat('_4d_exp_',tracer_type,'_gaussian_noise30');
img_type = 'synt'; % 'synt', 'it50', 'ground truth','real'
init_type = 'segmentation'; % type of initialization: 'nfindr','vca','ROI nfindr','ROI vca','rand','zero','ground truth','test A','test B','test M'
type_init = 'create'; % 'create','load'
K = 3;
V = 3;
if strcmp(tracer_type,'PE2I')
    tracer_decay=0.034; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
else if strcmp(tracer_type,'DPA')
%         tracer_decay=0.0063; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
        tracer_decay=0.06; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
    end
end
if strcmp(type_init,'create')
    
    % --------------------------------------------------------------
    % Choose image
    % --------------------------------------------------------------
    
    [PET,A_gt,B_gt,M_gt,alpha_gt,BP_gt,B_basic,vect_prop,M_ROI] = choose_image(img_type,img_name,K);
    
    keyboard
    %--------------------------------------------------------------
    % Initialization
    %--------------------------------------------------------------
    lbd = 0.5*max(var(PET.Y(:,PET.mask)));
    gamma= 0*max(var(PET.Y(:,PET.mask)));
    if strcmp(tracer_type,'PE2I')
        alpha_lim = [0.6,0.05];
        B_lim = [eps,-0.25,eps;0.7,-eps,0.15];
    else
        alpha_lim = [0.6,0.08];
        B_lim = [eps,-0.4,eps;0.7,-eps,0.15];
    end
    % % First estimation - initialization
    [A0,M0,B0,alpha0,BP0] = variables_initialization(init_type,PET,K,V,A_gt,lbd,gamma,alpha_lim,B_lim,tracer_decay);
    
    % % Put in the same order of GT for error computation
    % [A0,M0] = GT_compare(A0,M0,A_gt,M_gt,PET,img_type);
    
    save(strcat('mat\initializations\brain_init_',img_type,img_name))
    
else
    load(strcat('mat\initializations\brain_init_',img_type,img_name))
    % vect_prop = -vect_prop;
end

%--------------------------------------------------------------
% Unmixing parameters
%--------------------------------------------------------------

% PALM parameters
Niter = 3000;
beta = 0.5;   % endmember penalization parameter: 0.0086 - pen: 4.37, fct_term = 0.0916
lbd = 0.5; % variability penalization - sparsity
gamma= 0;
eta = 0.5; %abundance penalization - smoothing
if strcmp(tracer_type,'PE2I')
    alpha_lim = [0.6,0.05];
    B_lim = [eps,-0.25,eps;0.7,-eps,0.15];
else
    alpha_lim = [0.6,0.08];
    B_lim = [eps,-0.5,eps;0.7,-eps,0.15];
end
epsilon =5e-03;
C = max(var(PET.Y(:,PET.mask)));
%--------------------------------------------------------------
% Adjustments
%--------------------------------------------------------------

% A0 = A_gt;
% M0 = M_gt;
% B0 = B_gt;
% alpha0 = alpha_gt;

% A0(:,PET.mask) = rand(size(A_gt(:,PET.mask)));
% M0 = zeros(size(M_gt));
% B0(:,PET.mask,1) = eps*ones(size(B_gt(:,PET.mask,1)));
% B0(:,PET.mask,2) = -eps*ones(size(B_gt(:,PET.mask,2)));
% B0(:,PET.mask,3) = eps*ones(size(B_gt(:,PET.mask,3)));
% alpha0 = [0.6*ones(1,2);0.01*ones(1,2)];
%
% BP0 = B0(:,:,1);
% for k=1:K-1
%     BP0(k,:) = BP0(k,:)+sum(squeeze(B0(k,:,2:end))./repmat(alpha0(:,k)',size(A0,2),1),2)';
% end

%--------------------------------------------------------------
% DEPICT
%--------------------------------------------------------------
BP_depict = zeros(2,PET.W*PET.H*PET.D);
[BP_depict(1,:),alpha_depict,B_depict(:,:,1)] = depict(PET,M0(:,1),31,1e-3,lbd*max(var(PET.Y(:,PET.mask))),max(B_lim(:)),-max(B_lim(:)),0.03,0.6,500);
[BP_depict(2,:),alpha_depict,B_depict(:,:,2)] = depict(PET,M0(:,2),31,1e-3,lbd*max(var(PET.Y(:,PET.mask))),max(B_lim(:)),-max(B_lim(:)),0.03,0.6,500);
% %--------------------------------------------------------------
% % Algorithm
% %--------------------------------------------------------------
% A_palm(1,:) = zeros(size(A0(1,:)));
% A_palm(2:end,:) = A0;
% M_palm(:,1) = M_ROI;
% M_palm(:,2:end) = M0;
% B_palm = zeros(size(A1,:));
% [A_palm,M_palm,B_palm,elapsedTime_PALM,err,err_beta,err_A,err_M,obj_fct_palm,err_A1,ni_palm] = nonlin_PETunmixing(PET,A_palm,M_palm,A_gt,M_gt,B_gt,Niter,vect_prop,B_palm,eps,'PLMM',{M_ROI},'PENALTY M',{'DIFFERENCE',beta*C},'PENALTY A',{'SMOOTH',eta*C},'PENALTY B',{'SPARSE',lbd*C});


%--------------------------------------------------------------
% PNMM
%--------------------------------------------------------------

[A_cnu,M_cnu,B_cnu,alpha_cnu,BP_cnu,Q,elapsedTime,obj_fct,ni] = nonlin_PETunmixing(PET,A0,M0,B0,alpha0,Niter,lbd*C,gamma*C,beta*C,eta*C,alpha_lim,B_lim,epsilon,tracer_decay);
%--------------------------------------------------------------
% Compute SNR
%--------------------------------------------------------------
N_it = ni-1;
load(strcat('results\results_4D_GT_',tracer_type))
[err0,err0_MA] = cmpt_all_err_norm(PET,A0,M0,B0(:,:,1),alpha0,BP0,A_gt,M_gt,B_gt(:,:,1),alpha_gt,BP_gt,K,V);
[err_cnu,err_cnu_MA] = cmpt_all_err_norm(PET,A_cnu,M_cnu,B_cnu(:,:,1),alpha_cnu,BP_cnu,A_gt,M_gt,B_gt(:,:,1),alpha_gt,BP_gt,K,V);

% --------------------------------------------------------------
% Plot and save
%--------------------------------------------------------------


plot_PALM_brain_results_4D(PET,A_cnu,BP_cnu,M_cnu,A_gt,BP_gt,M_gt,A0,M0,BP0,BP_depict,obj_fct,K,img_type);
save(strcat('results\results_4D_',img_type,img_name))
% 1e-2
% err_cnu =
%     0.0489    0.0815    0.4226    0.0536
% 5e-3
% err_cnu =
%     0.0897    0.0858    0.4037    0.0803
