function [BP,alpha_dic,B] = onetissue_CM(PET,TAC_ref,V,epsilon,lbd,o_k2,u_k2,u_K1,o_K1,Niter)

% Dictionary of basis fuctions
k2_dic = log10(logspace(u_k2,o_k2));
V = length(K2_dic);

basis_fct = zeros(PET.L,V);
for Nd = 1:V
    E_toep =  toeplitz([exp(-PET.time'*k2_dic(Nd-1));zeros(PET.L-1,1)],[exp(-k2_dic(Nd-1)*PET.time(1));zeros(PET.L-1,1)]);
    basis_fct(:,Nd) = E_toep(1:PET.L,1:PET.L)*TAC_ref;
end
Y_tilde = PET.Y(:,PET.mask);

B = eps*ones(V,100);

for n=linspace(1,length(PET.mask),100)
    for i=1:V
        B(i,n) = basis_fct(:,i)\Y_tilde(:,n);
        allV(:,i,n) = basis_fct(:,i)*B(:,r,q);
        
    end
end

diff = squeeze(sum(abs(bsxfun(@minus,allV,Y_tilde)),1));

[minimum,ind] = min(diff(:));

K1_out = B(ind,n)
k2_out = k2_dic(ind)4