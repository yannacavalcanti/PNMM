function [ M] = palm_M_no_ROI(PET,A_old,M_old,V,B,M_ROI,M0,kernel3D,varargin)
% Compute endmembers estimation through PALM algorithm, except for ROI
% endmember
% 
%%%%%% Inputs:
% - Y               2D mixed image (L|N)
% - A_old           Initial abundances (K|N)
% - M_old           Initial Endmembers (L|K)
% - V               Internal eigenvectors (L|Nv)
% - B               Internal abundances (Nv|N)
% - beta           regularization parameter for M penalization term
% - M_ROI           mean of endmember ROI estimated from database
%%%%%% Outputs:
% - M               updated endmember
%
% Yanna Cavalcanti, Mars 2016


% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
% Exclude ROI endmember
M_est = M_old(:,2:end);
l = size(V,1);

% PSF convolution
AH =  psf_conv(A_old,PET.W,PET.H,PET.D,PET.mask,kernel3D);
delta = psf_conv(repmat(A_old(1,:),l,1).*(V*B),PET.W,PET.H,PET.D,PET.mask,kernel3D)+M_ROI*AH(1,:);

% Initialization
A = AH(2:end,PET.mask);
Y = PET.Y(:,PET.mask);
K = size(M_est,2);
Da = (delta(:,PET.mask)-Y)*A';
Ca = A*A';
Mp = zeros(K,K,K);
beta = varargin{1};

% Endmember penalization
%-----------------------------Mutual distance------------------------------
% for k = 1:K
%     Ek = zeros(K,1);
%     Ek(k) = 1;
%     Gk = -eye(K)+Ek*ones(1,K);
%     Mp(:,:,k) = Gk*Gk';
% end
% Mp = squeeze(sum(Mp,3));
% Reg = beta.*(M_est*Mp+2*(M_est-repmat(M_ROI,1,K)));
% Reg = beta.*(M_est*Mp);
%---------------------------L2 norm from start-----------------------------
Reg = beta*(M_est-M0(:,2:end));

mu = 1;

% L_endm = mu*(norm(Ca+Mp)+beta);
L_endm = mu*(norm(Ca)+beta);

grad_endm = M_est*Ca+Da+Reg;

% Compute update
M_new = M_est-(1/L_endm)*grad_endm;
% Project onto positive simplex
M_new(M_new<0) = 0;
M_new(M_new>1) = 1;

M = [M_ROI,M_new];