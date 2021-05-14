function [ M] = slmm_M(PET,A,M_old,V,B,M0,varargin)
% Compute endmember estimation through PALM algorithm
% 
%%%%%% Inputs:
% - Y               2D mixed image (L|N)
% - A               Initial abundances (K|N)
% - M_old           Initial Endmembers (L|K)
% - V               Internal eigenvectors (L|Nv)
% - B               Internal abundances (Nv|N)
% - beta           regularization parameter for M penalization term
%%%%%% Outputs:
% - M               updated endmember
%
% Yanna Cavalcanti, Mars 2016

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

l = size(V,1);
K = size(M_old,2);

% PSF convolution
AH =  A;
delta = repmat(A(1,:),l,1).*(V*B);

Y = PET.Y(:,PET.mask);
A = AH(:,PET.mask);
Da = (delta(:,PET.mask)-Y)*A';
Ca = A*A';
Mp = zeros(2,2,2);

beta = varargin{1};
% Endmember penalization
% for k = 1:2
%     Ek = zeros(2,1);
%     Ek(k) = 1;
%     Gk = -eye(2)+Ek*ones(1,2);
%     Mp(:,:,k) = Gk*Gk';
% end
% 
% Mp = beta.*squeeze(sum(Mp,3));
% 
% Mall=zeros(K,K);
% Mall(3:4,3:4)=Mp;
% Mp=Mall;

mu=1;
% Reg = M_old*Mp;
% L_endm = mu*(norm(Ca)+beta*norm(Mp));
Reg = beta*(M_old-M0);
Reg(:,1) = 100000*(M_old(:,1)-M0(:,1));
L_endm = mu*(norm(Ca)+beta);

grad_endm = M_old*Ca+Da+Reg;

% Compute update
M = M_old-(1/L_endm)*grad_endm;
M(:,1)=M_old(:,1);
% Project onto positive simplex
M(M<0) = 0;
