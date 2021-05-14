function [X] = four_to_threeD(X0,K,W,H,L)

X0 = reshape(X0,K,W,H,L);
X0 = X0(:,:,:,L/2);
X = reshape(X0,K,W*H);

