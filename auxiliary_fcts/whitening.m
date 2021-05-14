function [PET,W] = whitening(PET)    

I = PET.Y(:,PET.mask);

% % Estimate mean and covariance matrix
frame_mean = sum(I,2)/(length(PET.mask));
centeredData = bsxfun(@minus,I,frame_mean);
% covEst = centeredData*centeredData'/(length(PET.mask));
covEst = PET.cov;

% % Estimate variance - diagonal covariance matrix
% [c,l] = wavedec(PET.Y(:,PET.mask),1,'db3'); % wavelet decomposition
% noise_var = wnoisest(c,l,1:2); % find noise in the 2 first levels

% Standarization
[E,eigvalue,~] = svd(covEst);
D = sqrt(eigvalue);
W = D\E';
% W = D\ones(size(D));

% WHITEN THE frames
PET.Y(:,PET.mask) = W*centeredData+frame_mean;

%PET.Y(:,PET.mask) = W*centeredData+frame_mean;
% PET.Y(:,PET.mask) = W*PET.Y(:,PET.mask);

PET.Y(PET.Y<0)=0;