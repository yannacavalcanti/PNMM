
function PET = add_poisson_noise(PET)
PET.Y(:,PET.mask) = poissrnd(PET.Y(:,PET.mask)); %received signal
end