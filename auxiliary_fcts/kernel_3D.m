function G = kernel_3D(sigma_d, res)
    resol = res./min(res(:));
    siz = [ceil(sigma_d*3),ceil(sigma_d*3),ceil(sigma_d*3)];
    tama = ceil((siz-1)/2);
    if tama(1) < 1 | tama(2) < 1 |tama(3) < 1
        if sigma_d <= 0
            kernel = 1;
        else
            tama = [1 1 1];
            [x,y,z] = ndgrid(-tama(1):tama(1),-tama(2):tama(2),-tama(3):tama(3));

            kernel = exp(-(x.*x./(2.*(sigma_d./resol(1))^2) + ...
                y.*y/(2.*(sigma_d./resol(2))^2) + ...
                z.*z/(2.*(sigma_d./resol(3))^2)));
            G = kernel/sum(kernel(:));
        end
    else
        [x,y,z] = ndgrid(-tama(1):tama(1),-tama(2):tama(2),-tama(3):tama(3));
        kernel = exp(-(x.*x./(2.*(sigma_d./resol(1))^2) + ...
            y.*y/(2.*(sigma_d./resol(2))^2) + ...
            z.*z/(2.*(sigma_d./resol(3))^2)));
        G = kernel/sum(kernel(:));
    end
end