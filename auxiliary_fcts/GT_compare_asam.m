function [A0,M0] = GT_compare_asam(A0,M0,M)
K = size(A0,1);
err = zeros(K,K);
ord = zeros(1,K);
for i =1:K
    [err(i,:)] = asam(M(:,i),M0);
end
% Find minimum errors
[m,ind1] = min(err,[],1);
[m,ind2] = min(err,[],2);
% Find index in agreement
ind_agree=find(diag(ind1==ind2));
ord(ind_agree)=ind1(ind_agree);
ind1(ind_agree)=0;
ind2(ind_agree)=0;

ind1(sum(ind1==ord',1)~=0)=0;
ind2(sum(ind2==ord,2)~=0)=0;
ind_nonzero1=0;
ind_nonzero2=0;
% Find index which are nonzero for one and zero for other
while ~isempty(ind_nonzero1) || ~isempty(ind_nonzero2)
ind_nonzero1 = find(ind1~=0 & ind2'==0);
if ~isempty(ind_nonzero1)
ord(ind_nonzero1)=ind1(ind_nonzero1);
ind1(ind_nonzero1)=0;
end
ind_nonzero2 = find(ind2~=0 & ind1'==0);
if ~isempty(ind_nonzero2)
ord(ind_nonzero2)=ind2(ind_nonzero2);
ind1(ind_nonzero2)=0;
end
ind1(sum(ind1==ord',1)~=0)=0;
ind2(sum(ind2==ord,2)~=0)=0;
end
ind_nonzero = find(ind2~=0 & ind1'~=0);

for i=1:length(ind_nonzero)
    if ord(ind_nonzero(i))==0
        if ind1(ind_nonzero(i))~=0 & ind2(ind_nonzero(i))~=0
            ord(ind_nonzero(i))=min(ind1(ind_nonzero(i)),ind2(ind_nonzero(i)));
        else if ind1(ind_nonzero(i))~=0
                ord(ind_nonzero(i))=ind1(ind_nonzero(i));
            else if ind2(ind_nonzero(i))~=0
                    ord(ind_nonzero(i))=ind2(ind_nonzero(i));
                end
            end
        end
    end
ind1(sum(ind1==ord',1)~=0)=0;
ind2(sum(ind2==ord,2)~=0)=0;
end
%[m,indeq]=min([err(ind_nonzero,ind1(ind_nonzero)),err(ind2(ind_nonzero),ind_nonzero)]);
mis = setdiff(1:K,ord);
rep = find(sum(ord==ord')>1 | ord==0);
% 
% ord(min(rep))=mis;
abund=A0(ord,:);
endm = M0(:,ord);
A0=abund;
M0=endm;