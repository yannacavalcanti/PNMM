function [B,W] = pnmm_B(PET,A,M,B_old,W,Q,lbd,u_lim,o_lim)
% Compute abundance estimation through PALM algorithm
%
%%%%%% Inputs:
% - Y               2D mixed image (L|N)
% - A               Abundances (K|N)
% - M               Endmembers (L|K)
% - V               Internal eigenvectors (L|Nv)
% - B_old           Initial Internal abundances (Nv|N)
%%%%%% Outputs:
% - B               updated internal abundance
%
% Yanna Cavalcanti, Mars 2016

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

B = B_old;
mu = 1e-2; % regularization constant to give greater speed to the algorithm (may induce divergence)
A_tilde = A(1:end-1,PET.mask);
W = W(:,PET.mask,:);
V = size(B,3);

Y_tilde = PET.Y(:,PET.mask)-M*A(:,PET.mask);
for i=1:V
    Y_tilde = Y_tilde-Q(:,:,i)*W(:,:,i);
end

for i=1:V
    %Gradient
    grad_B = -(Q(:,:,i)'*Y_tilde).*A_tilde;
    %Lipschitz constant
    L_B = mu*norm(Q(:,:,i)'*Q(:,:,i))*norm(A_tilde.^2); % Norm majorisation
    
    % Compute update
    B_temp = B_old(:,PET.mask,i)-(1/L_B)*grad_B;
    
    % Group soft thr
    B_temp = proxL21col(B_temp,lbd/L_B)+eps;
    %Projection
    B_temp(B_temp<u_lim(i))=u_lim(i);
    B_temp(B_temp>o_lim(i))=o_lim(i);
    
    B(:,PET.mask,i) = B_temp;
end

W = repmat(A(1:end-1,:),1,1,V).*B;
