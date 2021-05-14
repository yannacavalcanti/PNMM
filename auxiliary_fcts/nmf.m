function [W,H] = nmf(V,W,H,epsilon,Niter)

err=1000;
i=2;
f = zeros(1,Niter);

while (err > epsilon) && (i<=Niter+1)

H = H.*(W'*V)./(W'*W*H);
W = W.*(V*H')./(W*(H*H'));

f(i) =  norm(V-W*H)^2;

err = abs(f(i-1)-f(i))/f(i-1);
i=i+1;

end