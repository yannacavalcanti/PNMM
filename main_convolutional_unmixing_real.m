clc,
clear all, format compact;
close all
addpath nonlinear_PETunmixing auxiliary_fcts;
addpath D:\ycruzcav\Database\mat
addpath mat plot\aux_fct plot\altmany-export_fig-5be2ca4

% --------------------------------------------------------------
% Choose initialization
% --------------------------------------------------------------
img_name = '_seg_02_02_3classes';
img_type = 'real'; % 'synt', 'it50', 'ground truth','real'
init_type = 'segmentation'; % type of initialization: 'nfindr','vca','ROI nfindr','ROI vca','rand','zero','ground truth','test A','test B','test M'
type_init = 'create'; % 'create','load'
version = '_v210512_test';
K = 3;
V = 3;
tracer_type = 'DPA';

if strcmp(tracer_type,'PE2I')
    tracer_decay=0.034; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
else if strcmp(tracer_type,'DPA')
%         tracer_decay=0.05; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
        tracer_decay=0.0063; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
    end
end
if strcmp(type_init,'create')
    % --------------------------------------------------------------
    % Load real image
    % --------------------------------------------------------------
    
    [PET,A0,M0] = load_real_image(K);
    
    %--------------------------------------------------------------
    % Initialization
    %--------------------------------------------------------------
    lbd = 0.2*max(var(PET.Y(:,PET.mask)));
    gamma= 0*max(var(PET.Y(:,PET.mask)));
    alpha_lim = [6,0.3];
    B_lim = [eps,-0.05,-0.01;1,0.3,0.05];

%     alpha_lim = [0.04,0.04];
%     B_lim = [-1,-1,-1;1,1,1];
    % % First estimation - initialization
    [A0,M0,B0,alpha0,BP0] = variables_initialization(init_type,PET,K,V,A0,lbd,gamma,alpha_lim,B_lim,tracer_decay);
    
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
beta = 0;   % endmember penalization parameter: 0.0086 - pen: 4.37, fct_term = 0.0916
lbd =2; % variability penalization - sparsity
gamma= 0;
eta = 0.7 ; %abundance penalization - smoothing
% eta=0.3; mu M=0.5 -- good version
% eta=0.8; mu M=0.45 -- very good version
alpha_lim = [6,0.3];
B_lim = [-1,-1,-1;3,3,3];
epsilon =5e-03;
% eta = 2.8; %abundance penalization - smoothing
% alpha_lim = [0.04,0.04];
% B_lim = [-1,-1,-1;1,1,1];
% eta = 0.15; %abundance penalization - smoothing
% epsilon =1e-03;
C = max(var(PET.Y(:,PET.mask)));

%--------------------------------------------------------------
% DEPICT
%--------------------------------------------------------------
BP_depict = zeros(2,PET.W*PET.H*PET.D);
[BP_depict(1,:),alpha_depict,B_depict(:,:,1)] = depict(PET,M0(:,1),31,epsilon,lbd*C,max(B_lim(:)),-max(B_lim(:)),tracer_decay,6,500);
[BP_depict(2,:),alpha_depict,B_depict(:,:,2)] = depict(PET,M0(:,2),31,epsilon,lbd*C,max(B_lim(:)),-max(B_lim(:)),tracer_decay,6,500);
%--------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------

[A_cnu,M_cnu,B_cnu,alpha_cnu,BP_cnu,Q,elapsedTime,obj_fct,ni] = nonlin_PETunmixing(PET,A0,M0,B0,alpha0,Niter,lbd*C,gamma*C,beta*C,eta*C,alpha_lim,B_lim,epsilon,tracer_decay);

% --------------------------------------------------------------
% Plot and save
%--------------------------------------------------------------

% plot_cnu_article_real_rotate(PET,A_cnu,BP_cnu,M_cnu,A0,M0,BP0,A0,M0,BP0,BP_depict,K,V)
save(strcat('results\results_4D_',img_type,img_name,version))
% 1e-2
% err_cnu =
%     0.0489    0.0815    0.4226    0.0536
% 5e-3
% err_cnu =
%     0.0897    0.0858    0.4037    0.0803
