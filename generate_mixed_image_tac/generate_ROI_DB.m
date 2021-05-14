function [vect_prop, B_mat,m_var] = generate_ROI_DB(L,Nv,A,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Database %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
vect_frame_time=[60 60 60 60 60 120 120 120 120 120 150 150 300 300 300 300 300 300 300 300]/64;
vect_frame_time = cumsum(vect_frame_time);

k3 = 0.001:1e-05:0.15;
% k3 = 0.001:0.01:0.15;

m_var = zeros(length(vect_frame_time),length(k3));

for i=1:length(k3)
    vect_parameters = [20 8 4 -1 -0.001 -0.001 0 0.25 0.4 k3(i) 0.1 0.01];
    % 
    % % parametres 5 et 6 sont ceux qui joue plus pour trouver la signature
    % % spécifique qui monte
    % % ja esta reglé

    m_var(:,i) = simulation_TAC_analytique(vect_parameters,vect_frame_time);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_var = m_var*max(M(:,1))/max(m_var(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = size(m_var,2);

Rmat = m_var-(mean(m_var,2)*ones(1,P));

Rmat = Rmat*Rmat'; % matrice de covariance empirique

[vect_prop D] = eigs(Rmat,Nv) ;
clear D
ind = find(A(1,:)~=0);
B_mat = repmat(A(1,:),L,1);
B_mat(:,ind) = m_var(:,unidrnd(size(m_var,2),1,length(ind)));
% B_mat = B_mat.*repmat(A(1,:),L,1);



