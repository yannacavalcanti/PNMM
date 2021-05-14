function [alpha,E,Q] = pnmm_alpha(PET,A,M,alpha_old,E,Q,W,gamma,t_toep,alpha_lim,tracer_decay)
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
A = A(:,PET.mask);
W = W(:,PET.mask,:);
alpha = alpha_old;
K = size(M,2);
V = size(alpha,1)+1;
mu = 4e-3; % regularization constant: <1 greater speed to the algorithm (may induce divergence)
% A tilde

Y0 = PET.Y(:,PET.mask)-M*A;
for i=1:V
    Y0 = Y0-Q(:,:,i)*W(:,:,i);
end

grad_alpha = zeros(V,K-1);
L_alpha = zeros(1,K-1);
for i=2:V
    G_alpha = zeros(K-1,K-1,K-1);
    for k=1:K-1
        Y_tilde = Y0+Q(:,k,i)*W(k,:,i);
        E_dev = (t_toep.*E(:,:,k,i-1));
        grad_alpha(i-1,k)= W(k,:,i)*(Y_tilde'*E_dev-(1/2)*W(k,:,i)'*M(:,k)'*(E_dev'*E(:,:,k,i-1)+E(:,:,k,i-1)'*E_dev))*M(:,k);
%         L_alpha(k) = norm(W(k,:,i))*(norm(-Y_tilde'+(1/2)*W(k,:,i)'*M(:,k)'*E(:,:,k,i-1)')+(3/2)*norm(W(k,:,i)'*M(:,k)')*norm(E(:,:,k,i-1)))*norm(E(:,:,k,i-1))*norm(t_toep)^2*norm(M(:,k));
        L_alpha(k) = norm(W(k,:,i))*(norm(-Y_tilde'+(1/2)*W(k,:,i)'*M(:,k)'*E(:,:,k,i-1)')+(3/2)*norm(W(k,:,i)'*M(:,k)')*norm(E(:,:,k,i-1)))*norm(E(:,:,k,i-1).*t_toep)*norm(t_toep)*norm(M(:,k));
%  -----------------------------Mutual distance----------------------------
        Ek = zeros(K-1,1);
        Ek(k) = 1;
        Gk = -eye(K-1)+Ek*ones(1,K-1);
        G_alpha(:,:,k) = Gk*Gk';
    end
    
    G_alpha = gamma*squeeze(sum(G_alpha,3));
    L_alpha = mu*(L_alpha+norm(G_alpha));
    grad_alpha(i-1,:) = grad_alpha(i-1,:)+alpha(i-1,:)*G_alpha;
    alpha(i-1,:) = alpha_old(i-1,:)-grad_alpha(i-1,:)./L_alpha;
    
    alpha(i-1,alpha(i-1,:)>alpha_lim(i-1))=alpha_lim(i-1);
end

alpha(alpha<tracer_decay)=tracer_decay;


for i=2:V
    for k=1:K-1
        E_toep =  toeplitz([exp(-PET.time'*alpha(i-1,k));zeros(PET.L-1,1)],[exp(-alpha(i-1,k)*PET.time(1));zeros(PET.L-1,1)]);
        E(:,:,k,i-1) = E_toep(1:PET.L,1:PET.L);
        Q(:,k,i) = E(:,:,k,i-1)*M(:,k);
    end
end