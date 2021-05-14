function [X] = proxL21col( U, gamma )
%PROXL21 proximal operator of L21 norm.
%  PROXL21  proximal operator of L21 norm using a image ordering
%
% Input:
% U:                optimization variable (3D image).
% gamma:            regularization parameter.

%
% Output:
% X:                 optimal point minimization.        
%
% AUTHOR:  Vinicius Ferraris
%          Institut de Recherche en Informatique de Toulouse
%          Departement Signal et Communications
%          Toulouse, 31000 France
%          vinicius.ferraris@enseeiht.fr
%
%  January 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gNorm = max(sqrt(sum(U.^2,1))-gamma,0); % L2 norm of columns
X = U.*repmat(gNorm./(gNorm+gamma),size(U,1),1);
