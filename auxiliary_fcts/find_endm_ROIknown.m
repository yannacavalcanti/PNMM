function [endm, matP, matU, Y_bar, endm_proj, y_proj] = find_endm_ROIknown(Y,R,endm_roi,method)


L_red = R-1; % number of principal component considered (endm. nbr - 1)
P = size(Y,2); %number of pixels
       
% PCA
disp('--> Begin - Principal Component analysis')
% [vect_prop y_red val_prop] = princomp(y);

Rmat = Y-(mean(Y,2)*ones(1,P));
Rmat = Rmat*Rmat'; % matrice de covariance empirique

OPTIONS.disp = 0; % diagnostic information display level
OPTIONS.maxit = 600; % maximum number of iterations

[vect_prop D] = eigs(Rmat,L_red,'LM',OPTIONS) ;
clear D
D = eye(L_red);
vect_prop = vect_prop';
disp('--> End - Principal Component analysis')

% first L_red eigenvectors
V = vect_prop(1:L_red,:);

% first L_red  eigenvalues
V_inv = pinv(V);
Y_bar = mean(Y,2);

% projector
matP =  D^(-1/2)*V; % permet de remédier aux problèmes numériques pouvant apparaître si le simplex est trop "étiré" 
                    % dans une direction (différence d'amplitude importante entre les valeurs propres)
% inverse projector
matU = V_inv*D^(1/2);

% projecting
y_proj = matP*(Y - Y_bar*ones(1,P));  % projection of the data on the R-1 principal axes.
t_roi = matP*(endm_roi-Y_bar); %projection of the ROI endmember


switch method
    case 'nfindr'
%            keyboard
        % NFINDR
        disp('--> Begin - N-FINDR')
        [endm_proj critere_opt] = nfindr_ROIknown(y_proj',t_roi);
        disp('--> End - N-FINDR')
          
        % in hyperspectral space
        endm = matU*endm_proj+Y_bar*ones(1,R); % inverse projection of the projected endm.

    case 'vca'
        
        % VCA
        disp('--> Begin - VCA')
        [endm] = vca_ROIknown(Y,endm_roi,R);
        disp('--> End - VCA')
        % in projected space
        endm_proj = matP*(endm -Y_bar*ones(1,R));
        
end