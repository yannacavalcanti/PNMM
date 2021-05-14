function [AHHt,AH,norm_HHt] = comp_smooth_grad(W,H,L,A,ord,Mask1)
% Compute smoothing norm with spatial penalization
% 
%%%%%% Inputs:
% - A           Initial abundances (K|N)
% - H,W     width and height of image
%%%%%% Outputs:
% - AHHt        Smoothing Penalization
% - norm_HHt    the norm of HHt
%
% Yanna Cavalcanti, April 2016

if ord==2 
    
[K,N] = size(A);
AHHt_h = zeros(K,N);
AHHt_v = zeros(K,N);
    % Compute H smoothing
    for i = 1:H
        AHHt_h(:,i) = A(:,i)-A(i+H);
        AHHt_h(:,i+((2:W-1)-1)*H) = 2*A(:,i+((2:W-1)-1)*H)-A(:,i+((2:W-1)-2)*H)-A(:,i+(2:W-1)*H);
        AHHt_h(:,(W-1)*H) = A(:,i+(W-1)*H)-A(:,i+(W-2)*H);
    end
    
    % Compute W smoothing
    for j = 1:W
        AHHt_v(:,1) = A(:,1)-A(:,2);
        AHHt_v(:,(2:H-1)+(j-1)*H) = 2*A(:,(2:H-1)+(j-1)*H)-A(:,(2:H-1)-1+(j-1)*H)-A(:,(2:H-1)+1+(j-1)*H);
        AHHt_v(:,(j)*H) = A(:,(j)*H)-A(:,(j)*H-1);
    end
    
    AHHt = AHHt_h+AHHt_v;
    AHHt = 2*AHHt;
    
    
    norm_HHt = 8; % (sqrt(2)+sqrt(2))^2
else
    

    [K,N] = size(A);
    A = reshape(A,K,W,H,L);
    I_Mask = zeros(size(A));
    I_Mask(:,Mask1) = 1;
    
    Mask = find(I_Mask==1);

    AHHt_h = NaN(size(A));
    AHHt_v = NaN(size(A));
    AHHt_d = NaN(size(A));
    
    AH_h = NaN(size(A));
    AH_v = NaN(size(A));
    AH_d = NaN(size(A));
    
    %%% Compute W smoothing
    % AH
    AH_v(:,1:end-1,:,:) = A(:,1:end-1,:,:)-A(:,2:end,:,:);
    ind = find(isnan(AH_v));
    AH_v(ind(ismember(ind,Mask))) = 0;
    % AHHt
    AHHt_v(:,1,:,:) = AH_v(:,1,:,:);
    AHHt_v(:,2:end-1,:,:) = AH_v(:,2:end-1,:,:) - AH_v(:,1:end-2,:,:);
    AHHt_v(:,end,:,:) = -AH_v(:,end-1,:,:);
    ind = find(isnan(AHHt_v));
%     AHHt_v(ind(ismember(ind,Mask))) = AH_v(ind(ismember(ind,Mask)));
    AHHt_v(ind(ismember(ind,Mask))) = 0;

    %%% Compute H smoothing
    % AH
    AH_h(:,:,1:end-1,:) = A(:,:,1:end-1,:)-A(:,:,2:end,:);
    ind = find(isnan(AH_h));
    AH_h(ind(ismember(ind,Mask))) = 0;

    % AHHt
    AHHt_h(:,:,1,:) = AH_h(:,:,1,:);
    AHHt_h(:,:,2:end-1,:) = AH_h(:,:,2:end-1,:) - AH_h(:,:,1:end-2,:);
    AHHt_h(:,:,end,:) = -AH_h(:,:,end-1,:);
    ind = find(isnan(AHHt_h));
%     AHHt_h(ind(ismember(ind,Mask))) = AH_h(ind(ismember(ind,Mask)));
    AHHt_h(ind(ismember(ind,Mask))) = 0;


    %%%% Compute L smoothing

    % AH
    if L>1
    AH_d(:,:,:,1:end-1) = A(:,:,:,1:end-1)-A(:,:,:,2:end);
    end
    ind = find(isnan(AH_d));
    AH_d(ind(ismember(ind,Mask))) = 0;

    % AHHt
    if L>1
    AHHt_d(:,:,:,1) = AH_d(:,:,:,1);
    AHHt_d(:,:,:,2:end-1) = AH_d(:,:,:,2:end-1) - AH_d(:,:,:,1:end-2);
    AHHt_d(:,:,:,end) = -AH_d(:,:,:,end-1);
    end
    ind = find(isnan(AHHt_d));
%     AHHt_d(ind(ismember(ind,Mask))) = AH_d(ind(ismember(ind,Mask)));
    AHHt_d(ind(ismember(ind,Mask))) = 0;

    %%% Total smoothing
    AHHt = AHHt_h+AHHt_v+AHHt_d;
    AH = AH_h+AH_v+AH_d;

    AH = reshape(AH,K,W*H*L);
    AHHt = reshape(AHHt,K,W*H*L);
    norm_HHt = 18; % (sqrt(2)+sqrt(2)+sqrt(2))^2 - spectral norm of 1st dim = sqrt(2) |[-1 1]x|/|x|
%     (N_dim*sqrt(2))^2
end
