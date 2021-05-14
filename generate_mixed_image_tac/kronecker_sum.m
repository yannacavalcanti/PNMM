function kro=kronecker_sum (m_label,i,j,K)

% function kro=kronecker_sum (m_label,i,j,K)
% Calculate the kronecker sum of the label prior
% 
% INPUTS
%      m_label  NcolxNrow             --- label matrix
%      i        1x1                   --- line indice
%      j        1x1                   --- colune indice
%      K        1x1                   --- number of classes
% 
% OUTPUT
%      kro      1x1                   --- result of the kronecker sum implemented


[Nrow, Ncol]= size(m_label);

kro=zeros(1,K);
% 4 neighboring pixels
v = [i+1 j; i-1 j; i j+1;i j-1];
% Check if pixels are in the first and last columns
    izero=[find(v(:,1)==0);find(v(:,2)==0);find(v(:,1)==Ncol+1);find(v(:,2)==Nrow+1)];
% Neglect values our of borders
    v(izero,:)=[];
%     Find value of neighboring chosen pixels
    z=diag(m_label(v(:,1),v(:,2)))';
%     Kronecker function - quantity of neighboring pixels belonging to each
%     class
    for ind=1:K
   kro(ind)=sum(z==ind);
    end



