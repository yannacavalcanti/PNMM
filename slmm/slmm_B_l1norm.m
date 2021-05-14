function [ B,L_beta] = slmm_B_l1norm(PET,A,M,V,B_old,lbd)
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

B = (B_old);
l = size(V,1);
Ea = repmat(A(1,:),l,1);
mu = 2; % regularization constant to give greater speed to the algorithm (may induce divergence)
D1 =M*A+Ea.*(V*B);
delta = (-PET.Y+D1);

% Computer gaussian filter norm
% H_norm = norm(fspecial('gaussian', [1 2*ceil(2*kernel3D)+1], kernel3D))^2;
% H_norm = 1;
%Lipschitz constant
L_beta = mu*(max(abs(A(1,:)))^2*norm(V)^2); % Norm majorisation

%Gradient
grad_beta = V'*(Ea(:,PET.mask).*delta(:,PET.mask));

% Compute update
B_new = B_old(:,PET.mask)-(1/L_beta)*grad_beta;

% Soft thr
B_new(abs(B_new)<=lbd/L_beta)=0;
B_new= sign(B_new).*(abs(B_new)-lbd/L_beta);

B_new(B_new<0)=0;
B(:,PET.mask)=B_new;

