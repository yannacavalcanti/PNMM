function [A0,M0] = GT_compare(A0,M0,A_gt,M_gt,PET,type)

K = size(A_gt,1);
ord = zeros(1,K);

abund = A0;
endm = M0;

if strcmp(type,'synt')
for i =1:K
    [diff,ind] = sort(sum(abs(M0-repmat(M_gt(:,i),1,K)),1));
    ord(i) = ind(1);
    rep=find(ord==ord(i));
    rep(rep==i)=[];
    if ~isempty(rep)
        rep(2)=i;
        [diff,ind] = sort(sum(abs(repmat(M0(:,ord(rep(1))),1,length(rep))-M_gt(:,rep)),1));
        mis = setdiff(1:K,ord);
        ord(rep(ind(2)))=mis(1);
    end

end
abund=A0(ord,:);
endm = M0(:,ord);

else
    [diff,ind] = min(sum(abs(A0(:,PET.SBR==1)-repmat(PET.SBR(PET.SBR==1),size(A0,1),1)),2));
%         [diff,ind] = min(sum(abs(M0-repmat(M_ROI,1,size(M0,2))),2));
abund(1,:) = A0(ind,:);
endm(:,1) = M0(:,ind);
abund(ind,:) = A0(1,:);
endm(:,ind) = M0(:,1);

end

A0 = abund;
M0 = endm;
