clc,
clear all, format compact;
close all


for num_img=1:20
    addpath palms auxiliary_fcts;
    addpath D:\ycruzcav\Database\mat
    addpath nonlinear_PETunmixing
    addpath mat
    % --------------------------------------------------------------
    % Choose initialization
    % --------------------------------------------------------------
    
    tracer_type = 'DPA';
    img_name = strcat('_4d_exp_',tracer_type,'_gaussian_noise20_',int2str(num_img));
    img_type = 'synt'; % 'synt', 'it50', 'ground truth','real'
    init_type = 'segmentation'; % type of initialization: 'nfindr','vca','ROI nfindr','ROI vca','rand','zero','ground truth','test A','test B','test M'
    type_init = 'load'; % 'create','load'
    K = 3;
    V = 3;
    if strcmp(tracer_type,'PE2I')
        tracer_decay=0.034; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
    else if strcmp(tracer_type,'DPA')
            tracer_decay=0.0063; %PE2I=0.034 (=0.693/half-life),DPA=0.0063
        end
    end
    if strcmp(type_init,'create')
        
        % --------------------------------------------------------------
        % Choose image
        % --------------------------------------------------------------
        
        [PET,A_gt,B_gt,M_gt,alpha_gt,BP_gt,B_basic,vect_prop,M_ROI] = choose_image(img_type,img_name,K);
        
        
        %--------------------------------------------------------------
        % Initialization
        %--------------------------------------------------------------
        lbd = 0.5*max(var(PET.Y(:,PET.mask)));
        gamma= 0*max(var(PET.Y(:,PET.mask)));
        if strcmp(tracer_type,'PE2I')
            alpha_lim = [6,0.05];
            B_lim = [eps,-0.25,eps;0.7,-eps,0.15];
        else
            alpha_lim = [6,0.01];
            B_lim = [eps,-0.4,eps;0.6,-eps,0.15];
        end
        % % First estimation - initialization
        [A0,M0,B0,alpha0,BP0] = variables_initialization(init_type,PET,K,V,A_gt,lbd,gamma,alpha_lim,B_lim,tracer_decay);
        
        % % Put in the same order of GT for error computation
        % [A0,M0] = GT_compare(A0,M0,A_gt,M_gt,PET,img_type);
        
        save(strcat('mat\initializations\brain_init_',img_type,img_name),'A0','M0','B0','alpha0','BP0')
        
    else
        % --------------------------------------------------------------
        % Choose image
        % --------------------------------------------------------------
        
        [PET,A_gt,B_gt,M_gt,alpha_gt,BP_gt,B_basic,vect_prop,M_ROI] = choose_image(img_type,img_name,K);
        
        
        
        load(strcat('mat\initializations\brain_init_',img_type,strcat('_4d_exp_',tracer_type,'_gaussian_noise30_',int2str(1))))
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
        alpha_lim = [6,0.05];
        B_lim = [eps,-0.25,eps;0.7,-eps,0.15];
    else
        alpha_lim = [6,0.08];
        B_lim = [eps,-0.5,eps;0.7,-eps,0.15];
    end
    epsilon =5e-03;
    C = max(var(PET.Y(:,PET.mask)));
    
    %--------------------------------------------------------------
    % DEPICT
    %--------------------------------------------------------------
    BP_depict = zeros(2,PET.W*PET.H*PET.D);
    [BP_depict(1,:),alpha_depict,B_depict(:,:,1)] = depict(PET,M0(:,1),31,1e-3,lbd*max(var(PET.Y(:,PET.mask))),max(B_lim(:)),-max(B_lim(:)),tracer_decay,6,500);
    [BP_depict(2,:),alpha_depict,B_depict(:,:,2)] = depict(PET,M0(:,2),31,1e-3,lbd*max(var(PET.Y(:,PET.mask))),max(B_lim(:)),-max(B_lim(:)),tracer_decay,6,500);
    
    % --------------------------------------------------------------
    % Plot and save
    %--------------------------------------------------------------
    
    
    % plot_PALM_brain_results_4D(PET,A_cnu,BP_cnu,M_cnu,A_gt,BP_gt,M_gt,A0,M0,BP0,BP_depict,obj_fct,K,img_type);
    save(strcat('results\DEPICT_results_',img_type,img_name),'BP_depict','alpha_depict','B_depict')
    clear all;
end