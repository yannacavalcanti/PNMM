clear all
close all
clc
addpath auxiliary_fcts;

err_all0 = zeros(7,20);
err_all_cnu = zeros(7,20);
err_all_palm = zeros(7,20);


mape_all0 = zeros(7,20);
mape_all_cnu = zeros(7,20);
mape_all_palm = zeros(7,20);

tracer_type = 'DPA';

for i=1:20
%     -------PALM
    load(strcat('results\PALM_results_4D_synt_4d_exp_',tracer_type,'_gaussian_noise20_',int2str(i)))
% NMSE
    err_all_palm(1:5,i)=err_palm;    
% MAPE
    
    [mape_palm(1)] = mape(A_palm(2:4,PET.mask),A_gt(:,PET.mask)); %A10
    [mape_palm(2)] = mape(M_palm(:,2:4),M_gt); %M0
    [mape_palm(3)] = 0;
    [mape_palm(4)] = 0;
    [mape_palm(5)] = 0;
    [mape_palm(6)] = 0;
    [mape_palm(7)] = 0;
    mape_all_palm(:,i)=mape_palm; 
%     -------PNMM
    load(strcat('results\results_4D_synt_4d_exp_',tracer_type,'_gaussian_noise20_',int2str(i)))
% NMSE
    err_all_cnu(1:5,i)=err_cnu;
    err_all_cnu(6,i) = cmpt_err_norm(BP_cnu(1,PET.mask),BP_gt(1,PET.mask));
    err_all_cnu(7,i) = cmpt_err_norm(BP_cnu(2,PET.mask),BP_gt(2,PET.mask));
% MAPE
    [mape_pnmm(1)] = mape(A_cnu(:,PET.mask),A_gt(:,PET.mask)); %A10
    [mape_pnmm(2)] = mape(M_cnu,M_gt); %M0
    [mape_pnmm(3)] = mape(B_cnu(:,PET.mask),B_gt(:,PET.mask));
    [mape_pnmm(4)] = mape(alpha_cnu,alpha_gt);
    [mape_pnmm(5)] = mape(BP_cnu(:,PET.mask),BP_gt(:,PET.mask));
    [mape_pnmm(6)] = mape(BP_cnu(1,PET.mask),BP_gt(1,PET.mask));
    [mape_pnmm(7)] = mape(BP_cnu(2,PET.mask),BP_gt(2,PET.mask));
    
    mape_all_cnu(:,i)=mape_pnmm;
%     -------Initialization
    load(strcat('mat\initializations\brain_init_synt_4d_exp_DPA_gaussian_noise20_',int2str(i)))
% NMSE    
    err_all0(1:5,i)=err0;
    err_all0(6,i) = cmpt_err_norm(BP0(1,PET.mask),BP_gt(1,PET.mask));
    err_all0(7,i) = cmpt_err_norm(BP0(2,PET.mask),BP_gt(2,PET.mask));
% MAPE
    [mape0(1)] = mape(A0(:,PET.mask),A_gt(:,PET.mask)); %A10
    [mape0(2)] = mape(M0,M_gt); %M0
    [mape0(3)] = mape(B0(:,PET.mask),B_gt(:,PET.mask));
    [mape0(4)] = mape(alpha0,alpha_gt);
    [mape0(5)] = mape(BP0(:,PET.mask),BP_gt(:,PET.mask));
    [mape0(6)] = mape(BP0(1,PET.mask),BP_gt(1,PET.mask));
    [mape0(7)] = mape(BP0(2,PET.mask),BP_gt(2,PET.mask));
    
    mape_all0(:,i)=mape0;
end
m0 = mean(err_all0,2);
m_cnu = mean(err_all_cnu,2);
m_palm = mean(err_all_palm,2);
all_m = [m0,m_cnu,m_palm]

v0 = var(err_all0,[],2);
v_cnu = var(err_all_cnu,[],2);
v_palm = var(err_all_palm,[],2);
all_v = [v0,v_cnu,v_palm]



mape_m0 = mean(mape_all0,2);
mape_m_cnu = mean(mape_all_cnu,2);
mape_m_palm = mean(mape_all_palm,2);
all_mape_m = [mape_m0,mape_m_cnu,mape_m_palm]


mape_v0 = var(mape_all0,[],2);
mape_v_cnu = var(mape_all_cnu,[],2);
mape_v_palm = var(mape_all_palm,[],2);
all_mape_v = [mape_v0,mape_v_cnu,mape_v_palm]
