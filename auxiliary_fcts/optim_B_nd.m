function [ B, B0, deltaM] = optim_B_nd(Y,A,M,ROI,vect_prop)

Nv = size(vect_prop,2);
% B = zeros(Nv,N);

E = Y - M*A;

W = vect_prop;
B0 = (W'*E);

B = B0./(ones(Nv,1)*A(ROI,:));

% deltaM = vect_prop*B;
deltaM=0;