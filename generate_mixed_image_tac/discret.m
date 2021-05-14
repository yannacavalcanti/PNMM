function Y = discret(K,p,N,M)

% function Y = discret(K,p,N,M)
% Execute the discretization of each pixel of the label vector according to
% the given probability, that is, put the more probable label in each pixel
% 
% INPUTS
%      K    1xNclass  --- vector with labels
%      p    1x1       --- probability
%      N    1x1       --- number of lines
%      M    1x1       --- number of colunes
% 
% OUTPUT
%      Y    1x1       --- label pixel resulted
%

%  Assure probability is normalized
p = p./sum(p);% p=vecteur de probabilite pour chaque classe

cum_p = [0 cumsum(p)]; %cummulative sum through classes
% cum_p(isnan(cum_p))=1;
NM = N*M;
z = rand(1,NM); % vecteur randomique de taille N*M

for i=1:NM
    ind = min(find(z(i)<cum_p)); % z(i)<cum_p: if z(i)< then the first probabilities, it means they are higher and the label must belong to that class
%   Define label for pixel
    Y(i) = K(ind-1);% on utilise la valeur de ind moins 1 parce qu'avec cum une dimension est rajouté
end
 
Y = reshape(Y,N,M); %reshape(A,m,n) returns the m-by-n matrix B whose 

