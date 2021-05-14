function [err_sep] = cmpt_err_sep(X,X0)
X =X';
X0 = X0';

D=abs(X'*X0);

err_sep = 20*log10((sum(sum(D./repmat(max(D,[],2),1,size(D,2)),2)-1)+sum(sum(D./repmat(max(D,[],1),size(D,1),1),1)-1))/(2*size(D,1)*(size(D,1)-1)));