function [Y,A0,M0,B,B_out,alpha,E,BP] = generate_parameters(A,M,V,vect_frame_time,K)
A0 = zeros(size(A));
[temp,ind] = max(A,[],1);
ind(sum(A)==0)=0;
% Clustering
for i=1:size(A0,1)
    A0(i,ind==i)=1;
end
A0(2,:) = A0(1,:)+A0(2,:);
A0(1,:) = [];
M0 = M(:,2:end);
% Specific binding in white matter
s = strel('sphere',4);
[s_x,s_y,s_z] = size(s.Neighborhood);
SBW = zeros(128,128,64);
SBW(72:72+s_x-1,75:75+s_y-1,27:27+s_z-1)=double(s.Neighborhood);
% Find maximum blood TAC of the image
[TACs_blood]=uniquetol((M*A(:,A0(3,:)==1))',0.5,'ByRows',true)';
M0(:,3) = TACs_blood(:,sum(TACs_blood)==max(sum(TACs_blood)));
% Find minimum gray matter TAC of the image
[TACs_GM]=uniquetol((M*A(:,A0(1,:)==1))',0.5,'ByRows',true)';
M0(:,1) = TACs_GM(:,sum(TACs_GM)==min(sum(TACs_GM)));

%Construct image
Y_tot = M*A;
ind_SBW=intersect(find((SBW==1)),find(A0(2,:)==1));
Y_tot(:,ind_SBW)=repmat(M(:,1),1,length(ind_SBW));
Y = M0*A0;

diff = sum(Y_tot-Y);
% Define parameters
vect_parameters_low = zeros(K-1,4);
vect_parameters_high = zeros(K-1,4);
vect_parameters_low(1,:) = [1,0.2,0.3,0.1];
vect_parameters_low(2,:) = [1,0.2,0.2,0.1];
vect_parameters_high(1,:) = [1.7,0.2,0.3,0.1];
vect_parameters_high(2,:) = [1.6,0.2,0.2,0.1];
% vect_parameters_high = [1.2,0.2,0.4,0.1];
% Initialize internal coefficients
B = zeros(K-1,length(diff),V);
B_out = zeros(K-1,2,V);
alpha = zeros(V-1,K-1);
for k=1:K-1
    % Put minimum TAC in all the indices smaller than it
    ind = intersect(find(diff<0.1),find(A0(k,:)==1));
    if ~isempty(ind)
        Y_tot(:,ind)=repmat(M0(:,k),1,length(ind));
    end
    % Define new TAC following the model for all the low uptake greater indices 
    ind1 = intersect(find(diff<150),find(diff>0.1)); 
    ind = intersect(ind1,find(A0(k,:)==1));
    [tac,alpha_out,B_out(k,1,:)] = FRTM_2tissue(vect_parameters_low(k,:),vect_frame_time,M0(:,k));
    if ~isempty(ind)
        alpha(:,k) = alpha_out;
        B(k,ind,:) = repmat(B_out(k,1,:),length(ind),1);
        Y_tot(:,ind)=repmat(tac,1,length(ind));
    end
        % Define new TAC following the model for all high uptake greater indices
    ind = intersect(find(diff>150),find(A0(k,:)==1));
    [tac,alpha_out,B_out(k,2,:)] = FRTM_2tissue(vect_parameters_high(k,:),vect_frame_time,M0(:,k));
    if ~isempty(ind)
%         alpha(:,k) = alpha_out;
        B(k,ind,:) = repmat(B_out(k,2,:),length(ind),1);
        Y_tot(:,ind)=repmat(tac,1,length(ind));
    end
end
k=K;
ind = intersect(find(diff<0.1),find(A0(k,:)==1));
if ~isempty(ind)
    Y_tot(:,ind)=repmat(M0(:,k),1,length(ind));
end


E = zeros(length(vect_frame_time),length(vect_frame_time),K-1,V-1);
Q = zeros(length(vect_frame_time),K-1,V);
Q(:,:,1) = M0(:,1:end-1);
for i=2:V
    for k=1:K-1
    E_toep =  toeplitz([exp(-vect_frame_time'*alpha(i-1,k));zeros(length(vect_frame_time)-1,1)],[exp(-alpha(i-1,k)*vect_frame_time(1));zeros(length(vect_frame_time)-1,1)]);
    E(:,:,k,i-1) = E_toep(1:length(vect_frame_time),1:length(vect_frame_time));
    Q(:,k,i) = E(:,:,k,i-1)*M0(:,k);
    end
end
Y = M0*A0;
for i=1:V
    Y = Y+Q(:,:,i)*(A0(1:end-1,:).*B(:,:,i));
end

BP = B(:,:,1);
for k=1:K-1
    BP(k,:) = BP(k,:)+sum(squeeze(B(k,:,2:end))./repmat(alpha(:,k)',length(diff),1),2)';
end
