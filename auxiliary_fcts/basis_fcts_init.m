function [B,alpha,BP,Q,elapsedTime,f,i] = basis_fcts_init(PET,A,M,B,alpha,Niter,lbd,gamma,alpha_lim,B_lim,epsilon,tracer_decay)
% pnmm unmixing and variables estimation
%
%%%%%% Inputs:
% - PET.Y           2D mixed image (PET.L|N)
% - PET.mask        mask of indices pointing to the real image
% - A               Abundances (K|N)
% - M               Endmembers (PET.L|K)
% - V               Internal eigenvectors (PET.L|Nv)
% - B               Internal abundances (Nv|N)
% - M_ROI           mean of endmember ROI estimated from database
% - PET.L               number of temporal slices
% - K               number of classes - 'nfindr','vca','ROI nfindr','ROI vca'
% - beta            regularization parameter for M penalization term
% - alpha           regularization parameter for A penalization term
% - PET.H,PET.W     width and height of image
%%%%%% Outputs:
% - A               updated abundance
% - M               updated endmember
% - B               updated internal abundance
% - elapsedTime     total processing time
%
% Yanna Cavalcanti, February 2016

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------

K = size(A,1);
M0 = M;
f = zeros(5,Niter+1);
V = size(B,3);
W = repmat(A(1:end-1,:),1,1,V).*B;
E = zeros(PET.L,PET.L,K-1,V-1);
Q = zeros(PET.L,K-1,V);
Q(:,:,1) = M(:,1:end-1);

for i=1:V-1
    for k=1:K-1
        E_toep =  toeplitz([exp(-PET.time'*alpha(i,k));zeros(PET.L-1,1)],[exp(-alpha(i,k)*PET.time(1));zeros(PET.L-1,1)]);
        E(:,:,k,i) = E_toep(1:PET.L,1:PET.L);
        Q(:,k,i+1) = E(:,:,k,i)*M(:,k);
    end
end
% Compute initial objective function
% f(:,1) = objective_pnmm(PET,M,A,B,alpha,Q,W,lbd,gamma,beta,eta,M0);
f(:,1) =1e10;
err = 1e10;

% Toeplitz matrix of time
t_toep = toeplitz([PET.time';zeros(PET.L-1,1)],[PET.time(1);zeros(PET.L-1,1)]);
t_toep = t_toep(1:PET.L,1:PET.L);
tic
i=2;
while (err > epsilon) & (i<=Niter+1)
    % Update B
    [B,W] = pnmm_B(PET,A,M,B,W,Q,lbd,B_lim(1,:),B_lim(2,:));

%     % Update alpha
    [alpha,E,Q] = pnmm_alpha(PET,A,M,alpha,E,Q,W,gamma,t_toep,alpha_lim,tracer_decay);

    % Objective function
    f(:,i)= objective_pnmm(PET,M,A,B,alpha,Q,W,lbd,gamma,0,0,M0);
    % Stopping criterion  or error
    err = (f(5,i-1)-f(5,i))/f(5,i-1)
    i =i+1
end

BP = B(:,:,1);
for k=1:K-1
    BP(k,:) = BP(k,:)+sum(squeeze(B(k,:,2:end))./repmat(alpha(:,k)',size(A,2),1),2)';
end

elapsedTime =  toc