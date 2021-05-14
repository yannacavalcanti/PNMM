function [ B] = slmm_B(PET,A,M,V,B_old)
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

B = NaN(size(B_old));
Y = PET.Y(:,PET.mask_nan);
A = A(:,PET.mask_nan);
l = size(V,1);
Ea = repmat(A(1,:),l,1);
eps = Y-M*A;
mu = 1; % regularization constant to give greater speed to the algorithm (may induce divergence)

% L_beta = ((norm(Ea))^2)*norm(V)^2;
% L_beta = mu*((max(abs(A(1,:)).^2))^2)*norm(V)^2; % Norm majorisation
L_beta = mu*(max(abs(A(1,:)))^2*norm(V)^2); % Norm majorisation

grad_beta = V'*(Ea.*(-eps+Ea.*(V*B_old(:,PET.mask_nan))));

% Compute update
B_new = B_old(:,PET.mask_nan)-(1/L_beta)*grad_beta;
B(:,PET.mask_nan)=B_new;