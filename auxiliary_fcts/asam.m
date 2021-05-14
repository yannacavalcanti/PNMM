function err = asam(X,X0)

% err = sum(acos(diag(X'*X0)/(norm(X)*norm(X0))))/size(X,2);
% err = sum(acos(sum(X.*X0,1)/(norm(X)*norm(X0))))/size(X,2);
K = min(size(X0));
err = zeros(1,K);
for i=1:K
err(i) = (acos(sum(X.*X0(:,i),1)/(norm(X)*norm(X0(:,i)))));
end