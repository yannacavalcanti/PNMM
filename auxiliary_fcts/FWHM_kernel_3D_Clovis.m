
function [G,G_norm] = FWHM_kernel_3D_Clovis(dim,type)

def_x=2.2; %mm
def_y=2.2; %mm
def_z=2.8; %mm

fwhm=4; %mm isotropic
% fwhm=6; %mm isotropic

sigma_g=fwhm/(2*sqrt(2*log(2))); % sigma gaussienne, en mm

%pour la taille du noyau on prend de façon à avoir plus de 6 sigma, en
%arrondissant au supérieur et pour un nombre de voxel impair
%de plus on veut un noyau homogene en dimension (de forme cubique), donc on
%se cale sur x ou y
taille=ceil(6*sigma_g/def_x);


noyau_x=zeros(taille,1);
noyau_y=zeros(taille,1);
noyau_z=zeros(taille,1);

%bornes inf et sup des intervalles pour le calcul des integrales
bornes_inf=zeros(taille,1);
bornes_sup=zeros(taille,1);

bornes_inf((taille+3)/2:end)=0.5:1:((taille-3)/2)+0.5;
bornes_sup((taille+1)/2:end)=0.5:1:((taille-1)/2)+0.5;

noyau_x=(erf(bornes_sup.*(def_x/sigma_g))-erf(bornes_inf.*(def_x/sigma_g)))./2;
noyau_y=(erf(bornes_sup.*(def_y/sigma_g))-erf(bornes_inf.*(def_y/sigma_g)))./2;
noyau_z=(erf(bornes_sup.*(def_z/sigma_g))-erf(bornes_inf.*(def_z/sigma_g)))./2;

noyau_x(1:(taille-1)/2)=noyau_x(end:-1:(taille+3)/2);
noyau_x((taille+1)/2)=2*noyau_x((taille+1)/2);
noyau_y(1:(taille-1)/2)=noyau_y(end:-1:(taille+3)/2);
noyau_y((taille+1)/2)=2*noyau_y((taille+1)/2);
noyau_z(1:(taille-1)/2)=noyau_z(end:-1:(taille+3)/2);
noyau_z((taille+1)/2)=2*noyau_z((taille+1)/2);
if dim == 3
    if strcmp(type,'3D')
        G=ones(taille,taille,taille);
        
        for i=1:taille
            for j=1:taille
                G(i,j,:)=squeeze(G(i,j,:)).*noyau_z;
            end
        end
        
        for j=1:taille
            for k=1:taille
                G(:,j,k)=squeeze(G(:,j,k)).*noyau_x;
            end
        end
        
        for i=1:taille
            for k=1:taille
                G(i,:,k)=squeeze(G(i,:,k)).*noyau_y';
            end
        end
    else
        G(1,:) = noyau_x;
        G(2,:) = noyau_y;
        G(3,:) = noyau_z;
    end
    G_norm = 2*norm(noyau_x)*norm(noyau_y)*norm(noyau_z);
    % G_norm = sqrt(sum(kernel3D(:).^2));
else
    %     G=ones(taille,taille);
    %
    %     for j=1:taille
    %         G(:,j)=squeeze(G(:,j)).*noyau_x;
    %     end
    %
    %     for i=1:taille
    %         G(i,:)=squeeze(G(i,:)).*noyau_y';
    %     end
    G(1,:) = noyau_x;
    G(2,:) = noyau_y;
    
    G_norm = norm(noyau_x)*norm(noyau_y);
    % G_norm = sqrt(sum(kernel3D(:).^2));
end