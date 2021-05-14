function AH = comp_smooth(A,H,W)
        [K,N] = size(A);
        
        %-Spatial smoothness constraint
        UpA = zeros(K,N);
        for k = 1:N-W
            UpA(:,W+k) = A(:,W+k) - A(:,k);
        end
        DownA  = [-UpA(:,W+1:end),zeros(K,W)];
        LeftA  = zeros(K,N);
        RightA = zeros(K,N);
        for k = 0:H-1
            LeftA(:,1+k*W:(k+1)*W)  = [zeros(K,1), diff(A(:,1+k*W:(k+1)*W),1,2)];
            RightA(:,1+k*W:(k+1)*W) = [-LeftA(:,2+k*W:(k+1)*W), zeros(K,1)];
        end
        AH = [LeftA,RightA,UpA,DownA];