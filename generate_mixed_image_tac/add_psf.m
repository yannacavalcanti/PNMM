function [PET] = add_psf(PET,type)
G = reshape(PET.Y,PET.L,PET.W,PET.H,PET.D);

switch type
    case '3D'
        def_x=2.2; %mm
        def_y=2.2; %mm
        def_z=2.8; %mm
        
        fwhm=4; %mm isotropic
        
        sigma_g=fwhm/(2*sqrt(2*log(2))); % sigma gaussienne, en mm
        sigma = [sigma_g/def_x,sigma_g/def_y,sigma_g/def_z];
        
        taille=ceil(6*sigma_g/def_x);
        for i=1:PET.L
            G(i,:,:,:) = imgaussfilt3(squeeze(G(i,:,:,:)),sigma,[taille,taille,taille]);
        end
        
    case '2D'
        
        for i=1:PET.L
            G(i,:,:,:) = imgaussfilt(squeeze(G(i,:,:,:)),sigma);
        end
        
    case '3D FWHM'
        
        [kernel3D,H_norm] = FWHM_kernel_3D_Clovis(3,'1D');
        for i=1:PET.L
            for j=1:PET.D
                %2D filt
                G(i,:,:,j) = conv2(kernel3D(1,:),kernel3D(2,:),squeeze(G(i,:,:,j)),'same');
            end
            G(i,:,:,:) = reshape(conv2(1,kernel3D(3,:),squeeze(reshape(G(i,:,:,:),PET.W*PET.H,PET.D)),'same'),PET.W,PET.H,PET.D);
        end
        
    case '2D FWHM'
        [kernel2D,H_norm] = FWHM_kernel_3D_Clovis(2,'1D');
        for i=1:PET.L
            %2D filt
            G(i,:,:) = conv2(kernel2D(1,:),kernel2D(2,:),squeeze(G(i,:,:)),'same');
        end
end

PET.Y(PET.Y>0) = G(PET.Y>0);