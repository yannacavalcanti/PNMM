clc,
clear all, format compact;
close all
addpath slmm auxiliary_fcts;
addpath D:\ycruzcav\Database\mat
addpath mat plot\aux_fct plot\altmany-export_fig-5be2ca4

% --------------------------------------------------------------
% Choose initialization
% --------------------------------------------------------------
img_name = '_seg_02_02_3classes';
img_type = 'real'; % 'synt', 'it50', 'ground truth','real'
init_type = 'segmentation'; % type of initialization: 'nfindr','vca','ROI nfindr','ROI vca','rand','zero','ground truth','test A','test B','test M'
type_init = 'create'; % 'create','load'
K = 3;
V = 3;
tracer_type = 'DPA';

if strcmp(tracer_type,'PE2I')
    tracer_decay=0.034; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
else if strcmp(tracer_type,'DPA')
        tracer_decay=0.0063; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
    end
end
if strcmp(type_init,'create')
    % --------------------------------------------------------------
    % Load real image
    % --------------------------------------------------------------
    
    [PET,A0,M0,M_ROI,vect_prop] = load_real_image(K);
    
    %--------------------------------------------------------------
    % Initialization
    %--------------------------------------------------------------
    lbd = 0.5*max(var(PET.Y(:,PET.mask)));
    gamma= 0*max(var(PET.Y(:,PET.mask)));
    alpha_lim = [0.6,0.11];
    B_lim = [-2,-2,-2;2,2,2];
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
beta = 0.1;   % endmember penalization parameter: 0.0086 - pen: 4.37, fct_term = 0.0916
lbd = 2.8; % variability penalization - sparsity
gamma= 0;
eta = 0.15; %abundance penalization - smoothing
alpha_lim = [0.05,0.01];
B_lim = [-2,-2,-2;2,2,2];
epsilon =1e-02;
C = max(var(PET.Y(:,PET.mask)));
%--------------------------------------------------------------
% Algorithm
%--------------------------------------------------------------
A_palm(1,:) = zeros(size(A0(1,:))); 
A_palm(2:4,:) = A0; 
M_palm(:,1) = M_ROI;
M_palm(:,2:4) = M0;
B_palm = zeros(size(A0(1,:)));
[A_palm,M_palm,B_palm,elapsedTime_PALM,err,err_beta,err_A,err_M,obj_fct_palm,err_A1,ni_palm] = slmm_unmix(PET,A_palm,M_palm,A0,M0,B0,Niter,vect_prop,B_palm,epsilon,'PLMM',{M_ROI},'PENALTY M',{'DIFFERENCE',beta*C},'PENALTY A',{'SMOOTH',eta*C},'PENALTY B',{'SPARSE',lbd*C});

% --------------------------------------------------------------
% Plot and save
%--------------------------------------------------------------

% plot_cnu_article_real_rotate(PET,A_palm,B_palm,M_palm,A0,M0,BP0,BP_depict,K,V)
save(strcat('results\PALM_results_4D_',img_type,img_name))
% 1e-2
% err_cnu =
%     0.0489    0.0815    0.4226    0.0536
% 5e-3
% err_cnu =
%     0.0897    0.0858    0.4037    0.0803
