function [ A] = slmm_A(PET,A_old,M,V,B,varargin)
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
l = size(V,1);
K = size(A_old,1);
W = V*B(:,PET.mask);
E = zeros(l,K);
E(:,1) = 1;
mu = 1; % regularization constant: <1 greater speed to the algorithm (may induce divergence)

% PSF Convolution
DH =M*A_old+repmat(A_old(1,:),l,1).*(V*B);

eps = (PET.Y-DH);

% Computer gaussian filter norm
% H_norm = norm(fspecial('gaussian', [1 2*ceil(2*kernel3D)+1], kernel3D))^2;
% H_norm = norm(kernel3D);

alpha = varargin{1};
% Matrix H
[AHHt,AH,norm_HHt] = comp_smooth_grad(PET.W,PET.H,PET.D,A_old,1,PET.mask);
clear AH
% norm_HHt = 0;
% AHHt = zeros(size(A_old));
%         [H,HHt,norm_HHt] = Nbr_operator(PET.W,PET.H,0,0);
%         clear H;
% L_abund = (norm(M'*M)+norm(E)*(2*norm(M)*(norm(W))+norm(W)^2)+alpha*l2Norm_DDt);

% Lipschitz constant
L_abund = mu*((norm(M'*M)+sqrt(l)*(2*norm(M)*(max(W(:))+max(W(:))^2))+alpha*norm_HHt)); % Norm majorisation
% Gradient
grad_abund = -M'*eps(:,PET.mask)-E'*(eps(:,PET.mask).*W)+alpha*AHHt(:,PET.mask);

% Compute update
A_new = A_old(:,PET.mask)-(1/L_abund)*grad_abund;

% %Sparsity
% [F,I] = sort(A_new,1,'descend');
% F(1:2,:) = proj_simplex_array(F(1:2,:));
% A_new = zeros(size(A_new));
% for k=1:2
%     for p=1:size(A_new,2)
%         A_new(I(k,p),p)=F(k,p);
%     end
% end
% A(:,PET.mask) = A_new;
% Project onto simplex

A(:,PET.mask) = proj_simplex_array(A_new);
