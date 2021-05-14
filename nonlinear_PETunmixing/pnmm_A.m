function [A,W] = pnmm_A(PET,A_old,M,Q,B_old,W,eta)
% Compute abundance estimation through PALM algorithm
% 
%%%%%% Inputs:
% - PET.Y           2D mixed image (L|N)
% - PET.mask_nan        mask_nan of indices pointing to the real image
% - A_old           Initial abundances (K|N)
% - M               Endmembers (L|K)
% - V               Internal eigenvectors (L|Nv)
% - B               Internal abundances (Nv|N)
% - L               number of temporal slices
% - K               number of classes - 'nfindr','vca','ROI nfindr','ROI vca'
% - alpha           regularization parameter for A penalization term
% - PET.H,PET.W     width and height of image
%%%%%% Outputs:
% - A               updated abundance
%
% Yanna Cavalcanti, Mars 2016

% -------------------------------------------------------------------------
% Project A into the simplex: sum(A,1) = 1 and A>0
% -------------------------------------------------------------------------
% The trick used here is that for every column of y, the result of bsxfun(...) is a vector whose elements are nondecreasing
% up to a certain index, at which the value is the desired threshold, and nonincreasing after. So, the desired threshold is simply the maximal element of this vector.
proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
A = A_old; 
K = size(A_old,1);
A_tilde = A_old(1:end-1,PET.mask);
M_tilde = M(:,1:end-1);
W = W(:,PET.mask,:);
B = B_old(:,PET.mask,:);
V = size(B,3);


mu = 5e-2; % regularization constant: <1 greater speed to the algorithm (may induce divergence)
% mu = 1e-2; % regularization constant: <1 greater speed to the algorithm (may induce divergence)

Y_tilde = PET.Y(:,PET.mask)-M*A_old(:,PET.mask);
for i=1:V
    Y_tilde = Y_tilde-Q(:,:,i)*W(:,:,i);
end

% Gradient
grad_A=-M'*Y_tilde;
L_A = norm(M'*M);
for i=1:V
    grad_A(1:end-1,:)=grad_A(1:end-1,:)-(Q(:,:,i)'*Y_tilde).*B(:,:,i);
    L_A = L_A+2*norm(M_tilde'*Q(:,:,i))*norm(B(:,:,i));
    for j=1:V
%         L_A = L_A+norm(B(:,:,i))*norm(Q(:,:,i)'*Q(:,:,j))*norm(B(:,:,j));
        L_A = L_A+norm(Q(:,:,i)'*Q(:,:,j))*norm(B(:,:,j).*B(:,:,i));
    end
end

% Matrix H
[AHHt,AH,norm_HHt] = comp_smooth_grad(PET.W,PET.H,PET.D,A_old,1,PET.mask);
clear AH
%Gradient
grad_A = grad_A+eta*AHHt(:,PET.mask);
% Lipschitz constant
L_A = mu*(L_A+eta*norm_HHt); % Norm majorisation

% Compute update
A_new = A_old(:,PET.mask)-(1/L_A)*grad_A;

A(:,PET.mask) = proj_simplex_array(A_new);

W = repmat(A(1:end-1,:),1,1,V).*B_old;
