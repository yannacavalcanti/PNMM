function [A,M,B,elapsedTime,err,err_B,err_A,err_M,f,err_A1,i] = slmm_unmix(PET,A,M0,A_gt,M_gt,B_gt,Niter,vect_prop,B,epsilon,varargin)
% slmm unmixing and variables estimation
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
f = zeros(5,Niter+1);
% err = zeros(1,Niter+1);
err = 1;
err_B = zeros(1,Niter+1);
err_A = zeros(1,Niter+1);
err_A1 = zeros(1,Niter+1);
err_M = zeros(1,Niter+1);

slmm_type = varargin{1};
M = M0;
for i=3:2:length(varargin)
    switch upper(varargin{i})
        case 'PENALTY M'
            varargin_objfct{i-2} = 'PENALTY M';
            type_M =  varargin{i+1}{1};
            switch upper(type_M)
                case 'MUTUAL DISTANCE'
                    beta = varargin{i+1}{2};
                    varargin_objfct{i-1} = {type_M,K,beta};
                case 'DIFFERENCE'
                    beta = varargin{i+1}{2};
                    varargin_objfct{i-1} = {type_M,M0,beta};
                case 'NONE'
                    beta=0;
                    varargin_objfct{i-1} = {type_M};
            end
            
        case 'PENALTY A'
            type_A = varargin{i+1}{1};
            varargin_objfct{i-2} = 'PENALTY A';
            switch upper(type_A)
                case 'SMOOTH'
                    alpha = varargin{i+1}{2};
                    varargin_objfct{i-1} = {type_A,alpha};
                case 'NONE'
                    alpha=0;
                    varargin_objfct{i-1} = {type_A};
            end
        case 'PENALTY B'
            type_B = varargin{i+1}{1};
            varargin_objfct{i-2} = 'PENALTY B';
            switch upper(type_B)
                case 'SPARSE'
                    lbd = varargin{i+1}{2};
                    varargin_objfct{i-1} = {type_B,lbd};
                    
                case 'NONE'
                    lbd = 0;
                    varargin_objfct{i-1} = {type_B};
            end
    end
end
var1 = {'PENALTY M',{'NONE'},'PENALTY A',{'NONE'},'PENALTY B',{'NONE'}};

% Compute initial objective function
f(:,1) = objective_p(PET,M,A,vect_prop,B,PET.L,var1);
f(:,1) = 1000;

tic
switch slmm_type
    case 'PLMM'
        M_ROI = varargin{2}{1};
        
        disp('--> Begin - ROI slmm')
        %         h = waitbar(0,'ROI slmm');
        i=2;
        while (err > epsilon) & (i<=Niter+1)

            % Update M
%             [ M] = slmm_M_no_ROI(PET,A,M,vect_prop,B,M_ROI,M0,kernel3D,beta);
            if i>50
            [M] = slmm_M(PET,A,M,vect_prop,B,M0,beta);
            end

            % Update A
            [A] = slmm_A(PET,A,M,vect_prop,B,alpha);
            
%             % Update B
            [ B] = slmm_B_l1norm(PET,A,M,vect_prop,B,lbd(1));
            
            f(:,i)= objective_p(PET,M,A,vect_prop,B,PET.L,varargin_objfct);

            % Stopping criterion  or error
            err = abs(f(5,i-1)-f(5,i))/f(5,i-1);
%             err = (f(5,i-1)-f(5,i))/f(5,i-1);
                i = i+1
%                 waitbar(i / (Niter+1))
            if i==51
                err=1;
            end
        end
%         % Compute B again - has to put the convolution
% 
%         B2=B;
%         B1=B;
%         ind = find(B(:,PET.mask)~=0);
%         B2(:,PET.mask(ind)) = (vect_prop'*(PET.Y(:,PET.mask(ind)) - M*A(:,PET.mask(ind))))./(ones(size(B,1),1)*A(1,PET.mask(ind)));
%         
%         for n = 1:length(ind)
%             
%             e = PET.Y(:,PET.mask(ind(n))) - M*A(:,PET.mask(ind(n)));
%             w = A(1,PET.mask(ind(n)))*vect_prop;
%             
%             B1(:,PET.mask(ind(n))) = ((w'*w)\(w'*e));
%             
%         end
%         B_nd = B;
%         [ B_nd(:,PET.mask(ind)), B0, deltaM] = optim_B_nd(PET.Y(:,PET.mask(ind)),A(:,PET.mask(ind)),M,1,vect_prop);
%              keyboard   
        err_A(i:end)  = [];
        err_A1(i:end) = [];
        err_M(i:end)  = [];
        err_B(i:end)  = [];
%         close(h)
        disp('--> End - ROI slmm')
        
    case 'STANDARD'
        disp('--> Begin - Standard slmm')
%         h = waitbar(0,'Standard slmm');
        i=2;
        while (err > epsilon) & (i<=Niter+1)
            
            % Update M
            [M] = slmm_M(PET,A,M,vect_prop,B,M0,kernel3D,beta);
            % Update A
            [ A] = slmm_A(PET,A,M,vect_prop,B,kernel3D,H_norm,alpha);
            
            f(:,i)= objective_p(PET,M,A,vect_prop,B,PET.L,kernel3D,lbdA,varargin_objfct);

            % Stopping criterion  or error
            err = abs(f(5,i-1)-f(5,i))/f(5,i-1);
                err_A(i)    = cmpt_err_norm(A(2:end,PET.mask),A_gt(2:end,PET.mask));
                err_A1(i)   = cmpt_err_norm(A(1,PET.mask),A_gt(1,PET.mask));
                err_M(i)    = cmpt_err_norm(M(:,2:end),M_gt(:,2:end));

               i = i+1;
%                 waitbar(i / (Niter+1))
        end
        err_A(i:end)  = [];
        err_A1(i:end) = [];
        err_M(i:end)  = [];
%         close(h)
        
        disp('--> End - Standard slmm')
end


elapsedTime =  toc