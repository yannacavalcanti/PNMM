function [ f ] = objective_p(PET,M,A,V,B,L,varargin)
% Compute objective function
% 
%%%%%% Inputs:
% - PET.Y           2D mixed image (L|N)
% - A               Abundances (K|N)
% - M               Endmembers (L|K)
% - deltaM          Variation V*B (L|N)
% - L               number of temporal slices
% - K               number of classes - 'nfindr','vca','ROI nfindr','ROI vca'
% - beta            regularization parameter for M penalization term
% - alpha           regularization parameter for A penalization term
% - PET.H,PET.w     width and height of image
%%%%%% Outputs:
% - f               objective function
%
% Yanna Cavalcanti, February 2016

varargin = varargin{1};

%--------------------------------------------------------------------------
% Initialization of terms
%--------------------------------------------------------------------------

for i=1:2:length(varargin)
    switch upper(varargin{i})
        case 'PENALTY M'
            %--------------------------------------------------------------
            % Endmember penalization term
            %--------------------------------------------------------------
            type =  varargin{i+1}{1};
            switch upper(type)
                case 'MUTUAL DISTANCE'
                    K = varargin{i+1}{2};
                    beta = varargin{i+1}{3};
                    for k = 1:K
                        Ek = zeros(K,1);
                        Ek(k) = 1;
                        Gk = -eye(K)+Ek*ones(1,K);
                        
                        Mp(:,:,k) = norm(M*Gk,'fro')^2;
                    end
                    termM = 0.5*beta*squeeze(sum(Mp,3));
                case 'DIFFERENCE'
                    M0 = varargin{i+1}{2};
                    beta = varargin{i+1}{3};
                    termM=0.5*beta*norm(M-M0)^2;
                case 'NONE'
                    termM=0;
            end
        case 'PENALTY A'
            %--------------------------------------------------------------
            % Abundance penalization term
            %--------------------------------------------------------------
            type = varargin{i+1}{1};
            switch upper(type)
                case 'SMOOTH'
                    alpha = varargin{i+1}{2};
                    A(isnan(A)) = 0;
                    [AHHt,AH,norm_HHt] = comp_smooth_grad(PET.W,PET.H,PET.D,A,1,PET.mask);
                    
                    clear AHHt norm_HHt
                    termA = 0.5*alpha*norm(AH(:,PET.mask),'fro')^2;
                case 'NONE'
                    termA = 0;
            end
        case 'PENALTY B'
            %--------------------------------------------------------------
            % Internal Abundance penalization term
            %--------------------------------------------------------------
            type = varargin{i+1}{1};
            switch upper(type)
                case 'SPARSE'
                    lbd = varargin{i+1}{2};
%                     termB = lbd(1)*norm(B(:,PET.mask),1)+0.5*lbd(2)*(B(:,PET.mask)*ones(length(PET.mask),1))^2;
                    termB = lbd(1)*norm(B(:,PET.mask),1);
                case 'NONE'
                    termB = 0;
            end
    end
end
%--------------------------------------------------------------
% Objective function
%--------------------------------------------------------------

% Objective

MAH = M*A+repmat(A(1,:),L,1).*(V*B);

f(1) = 0.5*(norm(PET.Y(:,PET.mask)-MAH(:,PET.mask),'fro')^2);
f(2) = termA;
f(3) = termM;
f(4) = termB;
f(5) = f(1)+termM+termA+termB;
