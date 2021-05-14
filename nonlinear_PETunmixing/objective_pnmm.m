function [ f ] = objective_pnmm(PET,M,A,B,alpha,Q,W,lbd,gamma,beta,eta,M0)
% Compute objective function
%
%%%%%% Inputs:
% - PET.Y           2D mixed image (L|N)
% - A               Abundances (K|N)
% - M               Endmembers (L|K)
% - deltaM          Variation V*B (L|N)
% - L               number of temporal slices
% - K               number of classes - 'nfindr','vca','ROI nfindr','ROI vca'
% - beta            regularization parameter for M penalization term
% - alpha           regularization parameter for A penalization term
% - PET.H,PET.w     width and height of image
%%%%%% Outputs:
% - f               objective function
%
% Yanna Cavalcanti, February 2016

V = size(B,3);
%--------------------------------------------------------------
% Endmember penalization term
%--------------------------------------------------------------

termM=0.5*beta*norm(M-M0)^2;

%--------------------------------------------------------------
% Abundance penalization term
%--------------------------------------------------------------

A(isnan(A)) = 0;
[AHHt,AH,norm_HHt] = comp_smooth_grad(PET.W,PET.H,PET.D,A,1,PET.mask);

clear AHHt norm_HHt
termA = 0.5*eta*norm(AH(:,PET.mask),'fro')^2;

%--------------------------------------------------------------
% Internal Abundance penalization term
%--------------------------------------------------------------
termB=0;
for i=1:V
termB = termB + lbd*sum(sqrt(sum(B(:,PET.mask,i).^2,1)));
end
%--------------------------------------------------------------
% Objective function
%--------------------------------------------------------------

% Objective

Y_tilde = M*A;
for i=1:V
    Y_tilde = Y_tilde+Q(:,:,i)*W(:,:,i);
end

f(1) = 0.5*(norm(PET.Y(:,PET.mask)-Y_tilde(:,PET.mask),'fro')^2);
f(2) = termA;
f(3) = termM;
f(4) = termB;
f(5) = f(1)+termM+termA+termB;
