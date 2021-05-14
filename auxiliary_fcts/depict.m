function [BP,alpha_dic,B] = depict(PET,TAC_ref,V,epsilon,lbd,o_lim,u_lim,u_alpha,o_alpha,Niter)

% V=31;
% Dictionary of basis fuctions
% alpha_dic = log10(logspace(u_alpha,o_alpha,V-1));
alpha_dic = log10(logspace(u_alpha,o_alpha,V));
alpha_dic = alpha_dic(2:end);

% TAC_ref = TAC_ref.*(1+exp(-PET.time'*u_alpha));

basis_fct = zeros(PET.L,V);
basis_fct(:,1) = TAC_ref;
for Nd = 2:V
        E_toep =  toeplitz([exp(-PET.time'*alpha_dic(Nd-1));zeros(PET.L-1,1)],[exp(-alpha_dic(Nd-1)*PET.time(1));zeros(PET.L-1,1)]);
        basis_fct(:,Nd) = E_toep(1:PET.L,1:PET.L)*TAC_ref;
end
Y_tilde = PET.Y(:,PET.mask)-repmat(TAC_ref,1,length(PET.mask));

B = eps*ones(V,length(PET.mask));
f(:,1) =1e10;
err = 1e10;
i=2;
while (err > epsilon) & (i<=Niter+1)
    grad_B = -basis_fct'*Y_tilde+basis_fct'*basis_fct*B;
    mu=1;
    L_B=mu*norm(basis_fct'*basis_fct);
    
    % Compute update
    B = B-grad_B/L_B;
    
    % Soft thr
    B = proxL21col(B,lbd/L_B);
    B = B+eps;
    
    %Projection
%     B(B>o_lim)=o_lim;
%     B(B<u_lim)=u_lim;
    
    % Objective function
    f(i)= 0.5*(norm(Y_tilde-basis_fct*B,'fro')^2)+lbd*sum(sqrt(sum(B.^2,1)),2);
    
    % Stopping criterion  or error
    err = (f(i-1)-f(i))/f(i-1);
    i =i+1
end

BP = nan(1,PET.W*PET.H*PET.D);
BP(:,PET.mask) = B(1,:)+sum(B(2:end,:)./repmat(alpha_dic',1,length(PET.mask)));