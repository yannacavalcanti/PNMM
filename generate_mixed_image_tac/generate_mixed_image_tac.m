function [PET,A,M,B,B_basic,alpha,E,BP] = generate_mixed_image_tac(PET,A0,M0,V,K,SNR,PSF,PSF_type,noise,tracer)

%%%%%%%%%%%%%%%%%%%%%%%%%% Generate image and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(tracer,'DPA')
    load data/dpa_tacs
    vect_frame_time = time_dpa;
    PET.L=length(vect_frame_time);
    [Y,A,M,B,B_basic,alpha,E,BP] = generate_parameters_dpa(A0,M0,V,vect_frame_time,K,M_dpa);
else
    % Number of f rames
    PET.L=size(M0,1);
    % Vector of time
    vect_frame_time=[60 60 60 60 60 120 120 120 120 120 150 150 300 300 300 300 300 300 300 300];
    vect_frame_time = cumsum(vect_frame_time)/60;
    [Y,A,M,B,B_basic,alpha,E,BP] = generate_parameters(A0,M0,V,vect_frame_time,K);
end
PET.Y = Y;

PET.mask = find(PET.Y(15,:)>0);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PSF==1
    [PET] = add_psf(PET,PSF_type);
end

% %%%%%%%%%%%%%%%%%%%%%%% Add white gaussian noise %%%%%%%%%%%%%%%%%%%%%%%%%%

if SNR>0
    switch noise
        case 'gaussian'
            PET = add_awgn_noise(PET,SNR,0);
        case 'gaussian_cov'
            PET = add_awgn_noise(PET,SNR,1);
        case 'poisson'
            PET = add_poisson_noise(PET);
        case 'gaussian_perpixel'
            PET.Y = PET.Y+sqrt(0.25*PET.Y).*rand(size(PET.Y));
    end
end

PET.time = vect_frame_time;