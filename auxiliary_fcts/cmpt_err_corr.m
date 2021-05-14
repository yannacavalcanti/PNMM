function [err_corr] = cmpt_err_corr(X,X0)
X=X';
X0=X0';

D=0;
for i=1:size(X,1)
D = D+(X(i,:)*X0(i,:)')/(norm(X(i,:))*norm(X0(i,:)));
end
err_corr=20*log10(-log10(D/size(X,1)));
