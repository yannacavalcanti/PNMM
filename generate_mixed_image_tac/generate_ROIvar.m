function [vect_prop, alpha, deltaM,m_var,M] = generate_ROIvar(L,Nv,SNR_ROI,A,M,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Database %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
vect_frame_time=[60 60 60 60 60 120 120 120 120 120 150 150 300 300 300 300 300 300 300 300]/64;
vect_frame_time = cumsum(vect_frame_time);

k3 = 0.008:1e-04:0.08;
% k3 = 0.001:1e-04:0.15;

m_var = zeros(length(vect_frame_time),length(k3));
if strcmp(type,'tac')
    for i=1:length(k3)
        vect_parameters = [20 8 4 -1 -0.001 -0.001 0 0.25 0.4 k3(i) 0.1 0.01];
        %
        % % parametres 5 et 6 sont ceux qui joue plus pour trouver la signature
        % % spécifique qui monte
        % % ja esta reglé
        
        m_var(:,i) = simulation_TAC_analytique(vect_parameters,vect_frame_time);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = size(m_var,2);
    
    Rmat = m_var-(mean(m_var,2)*ones(1,P));
    
    Rmat = Rmat*Rmat'; % matrice de covariance empirique
    
    [vect_prop D] = eigs(Rmat,Nv) ;
    clear D
    
    M_ROI = min(m_var,[],2);
%     M_ROI = (M_ROI-min(M_ROI))*(max(M(:,1))-min(M(:,1)))/(max(M_ROI)-min(M_ROI))+min(M(:,1));
    m_var = (m_var-min(M_ROI))*(max(M(:,1))-min(M(:,1)))/(max(M_ROI)-min(M_ROI))+min(M(:,1));

    M(:,1) = min(m_var,[],2);
    % clear D

else if strcmp(type,'gauss') vect_prop = normpdf(1:20,12,3)'; M_ROI = M(:,1);
    else
        vect_prop = [zeros(10,1);ones(10,1)]; M_ROI = M(:,1);
    end
end

vect_prop = abs(vect_prop);
% 
% 
% % Random B
% B_mat = A(1,:);
% % ind = find(~isnan(B_mat));
% 
% B_mat(B_mat<0.5) = 0;
% 
% ind = find(B_mat>0.5);
% 
% B_mat = B_mat(:);
% 
% % Signal Power
% SNR = 10^(SNR_ROI/10); %SNR to linear scale
% P_m = sum(abs(M_ROI).^2)/(L); %Calculate actual symbol energy
% P_v = sum(abs(vect_prop(:)).^2)/(Nv*L);
% P_b = P_m/P_v;
% C_b=P_b/SNR; %Find the noise spectral density
% 
% alpha = repmat(B_mat',Nv,1);
% alpha(:,ind) = sqrt(C_b)*randn(Nv,length(ind));
% deltaM = vect_prop(:,1:Nv)*alpha;
% 
% 
% 

% B with regions
B_mat = reshape(A(1,:),128,128,64);

B_mat(B_mat<0.5) = 0;

SNR_ROI = 12; % to match the minimum
% Signal Power
SNR = 10^(SNR_ROI/10); %SNR to linear scale
P_m = sum(abs(M(:,1)).^2)/(L); %Calculate actual symbol energy
P_v = sum(abs(vect_prop(:)).^2)/(Nv*L);
P_b = P_m/P_v;
C_b=P_b/SNR; %Find the noise spectral density

alpha = repmat(B_mat,Nv,1);
alpha1 = repmat(B_mat(:,1:64,:),Nv,1);
alpha2 = repmat(B_mat(:,64:end,:),Nv,1);
alpha1(alpha1>0.9 & alpha1<0.95) = sqrt(C_b)/2+sqrt((sqrt(C_b)/2)/SNR)*randn(size(find(alpha1>0.9 & alpha1<0.95)));
alpha1(alpha1>0.95 & alpha1<1) = sqrt(C_b)+sqrt((sqrt(C_b))/SNR)*randn(size(find(alpha1>0.95 & alpha1<1)));
alpha2(alpha2>0.9 & alpha2<0.95) = 3*sqrt(C_b)/2+sqrt((3*sqrt(C_b)/2)/SNR)*randn(size(find(alpha2>0.9 & alpha2<0.95)));
alpha2(alpha2>0.95 & alpha2<1) = 2*sqrt(C_b)+sqrt(2*sqrt(C_b)/SNR)*randn(size(find(alpha2>0.95 & alpha2<1)));
alpha(:,1:64,:) = alpha1;
alpha(:,64:end,:) = alpha2;

alpha = reshape(alpha,1,128*128*64);
deltaM = vect_prop(:,1:Nv)*alpha;
