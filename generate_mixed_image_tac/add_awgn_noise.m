 
function [PET] = add_awgn_noise(PET,SNR_dB,factors)
x = PET.Y(:,PET.mask);
%y=awgn_noise(x,SNR) adds AWGN noise vector to signal 'x' to generate a
%resulting signal vector y of specified SNR in dB
rng('default');%set the random generator seed to default (for comparison only)
L=length(x(:));
SNR = 10^(SNR_dB/10); %SNR to linear scale
Esym=sum(abs(x(:)).^2)/(L); %Calculate actual symbol energy
N0=Esym/SNR; %Find the noise spectral density
%noiseSigma = sqrt(N0);%Standard deviation for AWGN Noise when x is real
noiseSigma = N0;%
n = randn(size(x));%computed noise

if factors==1
    load D:\ycruzcav\Database\mat\noise_factors
%     noise_factors = normpdf(1:20,1,3)';
    noise_factors = (noise_factors*0.9/max(noise_factors)+0.1);
%     noise_factors = (noise_factors*0.9/max(noise_factors)+0.1).*rand(20,1);
    PET.cov=noiseSigma*eye(size(x,1),size(x,1)).*repmat(noise_factors,1,length(noise_factors));

    [V,D] = eig(PET.cov);
    n = (V*sqrt(D))*n;
else
    n = sqrt(noiseSigma)*n;
end
 PET.Y(:,PET.mask) = x + n; %received signal
% PET.Y(PET.Y<0)=0;