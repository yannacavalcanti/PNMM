function [endm,critere_opt] = nfindr_ROIknown(y_proj,t_roi)

% Modified by Yanna Cavalcanti: Jan 2016


% knowing y_proj = D^(-1/2)*V*(Y - Y_bar*ones(1,P)) and U = V_inv*eye(L_red)^(1/2)
% size(y_red) = Nb_pix x Nb_bandes
[Nb_pix,Nb_bandes] = size(y_proj);


% Find the difference between the projection of ROI endmember and the
% signature of the pixels so we can delete the equivalent pixels. lim is
% the difference allowed to consider they are equivalent
lim = 0.3;
t1 = find(mean(abs(y_proj-repmat(t_roi',Nb_pix,1)),2)<lim);
[G,ind] = sort(sum(abs(y_proj(t1,:)-repmat(t_roi',length(t1),1)),2));
clear G
t2 = [];
if ~isempty(ind)
t2 = t1(ind(1)); % index of the projection correspondent to the ROI endm
end
y_proj(t2,:) = [];

[Nb_pix,Nb_bandes] = size(y_proj);

Nb_endm = Nb_bandes+1;

% quick hull if size allows it
if Nb_endm<7
%   K = CONVHULLN(X) returns the indices K of the points in X that 
%   comprise the facets of the convex hull of X. 
    tmp = convhulln(y_proj);
%   C = UNIQUE(A) for the array A returns the same values as in A but with 
%   no repetitions. C will be sorted.  
    ind_envlp = unique(tmp(:));
else
    ind_envlp = 1:Nb_pix;
end

envlp = y_proj(ind_envlp,:);
Nb_pix_envlp = length(ind_envlp);

% INITIALIZATION
% choice of Nb_endm pixels
%   P = RANDPERM(N) returns a vector containing a random permutation of the
%   integers 1:N.  For example, RANDPERM(6) might be [2 4 5 6 1 3].
ind_perm = randperm(Nb_pix_envlp);
combi = sort(ind_perm(1:Nb_endm-1));

% cacul du critère

candidat_opt = [t_roi,y_proj(combi,1:Nb_bandes)']; %add endm proj as a cadidate
% volume calculation - the denominator is always the same
critere_opt = -abs(det([candidat_opt;ones(1,Nb_endm)]));

% Algorithm
h = waitbar(0,['N-FINDR...'],'CreateCancelBtn','closereq', 'name', 'N-FINDR');
for ind_endm=1:Nb_endm-1
    for ind_pix=1:Nb_pix_envlp
        waitbar((Nb_pix_envlp*(ind_endm-1)+ind_pix)/(Nb_pix_envlp*Nb_endm),h)

        combi_cand = combi;
%         Replacement: for each pixel vector, recalculate the volume by
%         testing the pixel in all p endmember positions
        combi_cand(ind_endm) = ind_pix;

        if length(unique(combi_cand)) == Nb_endm-1 % i~=j

            %             Envelop candidate
            candidat_cand = [t_roi,envlp(combi_cand,1:Nb_bandes)']; %add endm proj as a cadidate

%             Volume calculation
            critere_cand = -abs(det([candidat_cand;ones(1,Nb_endm)]));

            %test
            if critere_cand<critere_opt %because it is negative
                critere_opt = critere_cand;
                candidat_opt = candidat_cand;
                combi = combi_cand;
            end
        end
    end
end
delete(h)
endm = candidat_opt;



