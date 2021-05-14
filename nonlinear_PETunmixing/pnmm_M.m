function [ M,Q] = pnmm_M(PET,A,M_old,Q,E_old,W,M0,beta)
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
M = M_old;
K = size(M,2);
W = W(:,PET.mask,:);
A = A(:,PET.mask);
V = size(W,3);
E = zeros(PET.L,PET.L,K-1,V);
E(:,:,:,1) = repmat(eye(PET.L),1,1,K-1);
E(:,:,:,2:end) = E_old;

Y0 = PET.Y(:,PET.mask)-M*A;
for i=1:V
    Y0 = Y0-Q(:,:,i)*W(:,:,i);
end

grad_M = zeros(size(M));
% L_M = zeros(1,K);
L_M= norm(A*A')+beta;
for k=1:K-1
    Y_tilde = Y0+M(:,k)*A(k,:);
    for i=1:V
        Y_tilde = Y_tilde+Q(:,k,i)*W(k,:,i);
    end
    grad_M(:,k) = -(Y_tilde-M(:,k)*A(k,:))*A(k,:)';

    for i=1:V
        Eki = E(:,:,k,i);
        grad_M(:,k) = grad_M(:,k)-Eki'*Y_tilde*W(k,:,i)'+(Eki+Eki')*M(:,k)*W(k,:,i)*A(k,:)';
        L_M = L_M+norm(Eki+Eki')*norm(W(k,:,i)*A(k,:)');
        for j=1:V
            Ekj = E(:,:,k,j);
            grad_M(:,k) = grad_M(:,k)+(1/2)*(Eki'*Ekj+(Eki'*Ekj)')*(M(:,k)*W(k,:,j)*W(k,:,i)');
            L_M = L_M+norm(Eki'*Ekj)*norm(W(k,:,j)*W(k,:,i)');
        end
    end
end

grad_M(:,K) = -Y0*A(K,:)';
% L_M(K) = norm(A(:,K)*A(:,K)')+beta;

% mu=0.5; % official
mu=0.45; %test
Reg = beta*(M_old-M0);

% Gradient
grad_M = grad_M +Reg;
%Lipschitz constant
L_M = mu*L_M;

% Compute update
% M = M_old-grad_M./repmat(L_M,PET.L,1);
M = M_old-grad_M./L_M;

% Project onto positive simplex
M(M<0) = 0;

Q = zeros(PET.L,K-1,V);
Q(:,:,1) = M(:,1:end-1);
for i=2:V
    for k=1:K-1
        Q(:,k,i) = E_old(:,:,k,i-1)*M(:,k);
    end
end
