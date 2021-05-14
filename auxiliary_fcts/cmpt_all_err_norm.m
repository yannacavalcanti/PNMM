function [err0,err0_MA] = cmpt_all_err_norm(PET,A0,M0,B0,alpha0,BP0,A_gt,M_gt,B_gt,alpha_gt,BP_gt,K,V)
% B0 = reshape(B0(:,PET.mask,:),K-1,length(PET.mask)*V);
% B_gt = reshape(B_gt(:,PET.mask,:),K-1,length(PET.mask)*V);
%%%%%%%%%%%%%%%%%%% Sunsal/Nfindr
[err0(1)] = cmpt_err_norm(A0(:,PET.mask),A_gt(:,PET.mask)); %A10
[err0(2)] = cmpt_err_norm(M0,M_gt); %M0
[err0(3)] = cmpt_err_norm(B0(:,PET.mask),B_gt(:,PET.mask));
[err0(4)] = cmpt_err_norm(alpha0,alpha_gt);
[err0(5)] = cmpt_err_norm(BP0(:,PET.mask),BP_gt(:,PET.mask));

% snr_M10 = sum(sqrt(sum(abs(M1_hat0(:,PET.mask)-B_mat(:,PET.mask)).^2,1)))/length(PET.mask);

% MA
for i=1:size(A0,1)
    [err0_MA(i)] = cmpt_err_norm(M0(:,i)*A0(2,PET.mask),M_gt(:,i)*A_gt(2,PET.mask));
end
