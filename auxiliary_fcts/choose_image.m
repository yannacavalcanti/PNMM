function [PET,A_gt,B_gt,M_gt,alpha_gt,BP_gt,B_basic,vect_prop,M_ROI] = choose_image(type,type1,K)

% Choose the image to be estimated - synthetic, ground truth, ...

switch type
    case 'real'
        load('D:\ycruzcav\tours\PET_unmixing\mat\initializations\brain_init_real_seg_02_03_3classes','A0','PET','M0')
        abund = reshape(permute(reshape(A0(1:2,:)',PET.W,PET.H,PET.D,K-1),[2 1 3 4]),PET.W*PET.H*PET.D,K-1)';
        endm = M0(:,1:2);

        load(strcat('D:\ycruzcav\ycruzcav_git\pet_unmixing\code\PET_unmixing_psf\mat\initializations\brain_init_real_5class'),'A0','A_gt','B0','B_gt','M0','M_gt','PET')
        abund2 = A0;
        endm2 = M0;
        mask2 = find(A0(5,:)==1);
        
        z = zeros(1,size(PET.Y,2));
        z(PET.mask)=1;
        z(mask2)=0;
        PET.mask = find(z==1);
        K=3;
        clear A0 M0
        A0(1:2,:)=abund;
        A0(3,:)=abund2(3,:);
        
        M0(:,1:2)=endm;
        M0(:,3)=endm2(:,3);
        vect_frame_time=[10*ones(1,6) 30*ones(1,8) 60*ones(1,4) 120*ones(1,5) 300*ones(1,8)];
        PET.time = cumsum(vect_frame_time);
        alpha_gt = 0;
        BP_gt = 0;
        B_basic = 0;
    case 'varPET'
        load pet_img_sample_var
        
        %         M_ROI= min(m_var,[],2);
        PET.C = max(PET.Y(:));
        PET.Y = PET.Y/PET.C;
        M_ROI = M_ROI/PET.C;
        vect_prop = vect_prop(:,1);
        Nv = size(vect_prop,2);
        A_gt = zeros(K,size(PET.Y,2));
        M_gt = zeros(PET.L,K);
        B_gt = zeros(Nv,size(PET.Y,2));
        B_mat = zeros(size(PET.Y));
        vect_prop = vect_prop(:,1);
    case 'it06'
        load pet_img_sample_it06
        load mat\brain_mask
        PET.mask  = PET.mask;
        
        PET.C = max(PET.Y(:));
        PET.Y = PET.Y/PET.C;
        M_ROI = M_ROI/PET.C;
        vect_prop = vect_prop(:,1);
        Nv = size(vect_prop,2);
        A_gt = zeros(K,size(PET.Y,2));
        M_gt = zeros(PET.L,K);
        B_gt = zeros(Nv,size(PET.Y,2));
        B_mat = zeros(size(PET.Y));
        vect_prop = vect_prop(:,1);
    case 'it16'
        load pet_img_sample_it16
        load mat\brain_mask
        PET.mask  = PET.mask;
        
        PET.C = max(PET.Y(:));
        PET.Y = PET.Y/PET.C;
        M_ROI = M_ROI/PET.C;
        vect_prop = vect_prop(:,1);
        Nv = size(vect_prop,2);
        M_gt = zeros(PET.L,K);
        B_gt = zeros(Nv,size(PET.Y,2));
        B_mat = zeros(size(PET.Y));
        vect_prop = vect_prop(:,1);
        PET.SBR = A_gt(1,:);
        PET.SBR(PET.SBR>0.3)=1;
    case 'it50'
        load pet_img_sample_50it
        PET.Y = PET.Y;
        i = find(PET.Y(15,:)>3.5e04);
        PET.mask = find(sum(PET.Y,1)>1e5);
        PET.H = H;
        PET.W = W;
        PET.D=1;
        
        %--------------------------------------------------------------
        % Find mean ROI
        %--------------------------------------------------------------
        m_var = PET.Y(:,i);
        
        M_ROI = mean(m_var,2);
        PET.C = max(PET.Y(:));
        PET.Y = PET.Y/PET.C;
        M_ROI = M_ROI/PET.C;
        vect_prop = vect_prop(:,1);
    case 'ground truth'
        load pet_ground_truth
        PET.Y = PET.Y;
        i = find(PET.Y(15,:)>36.05);
        PET.mask = find(sum(PET.Y,1)>50);
        PET.H = H;
        PET.W = W;
        PET.D=1;
        
        %--------------------------------------------------------------
        % Find mean ROI
        %--------------------------------------------------------------
        m_var = PET.Y(:,i);
        
        M_ROI = mean(m_var,2);
        PET.C = max(PET.Y(:));
        PET.Y = PET.Y/PET.C;
        M_ROI = M_ROI/PET.C;
        vect_prop = vect_prop(:,1);
    case 'synt'
        im_name = strcat('mat\generated_images\brain_synt',type1);
        load(im_name)
        % Set NaNs
        Y_NaN = NaN(size(PET.Y));
        Y_NaN(:,PET.mask) = PET.Y(:,PET.mask);
        A_gt_NaN = NaN(size(A_gt));
        A_gt_NaN(:,PET.mask) =A_gt(:,PET.mask);
        A_gt = A_gt_NaN;
        clear A_gt_NaN;
        PET.Y = Y_NaN;
        % Eigenvector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind = find(sum(BP_gt)~=0 & ~isnan(sum(BP_gt)));
        m_var = PET.Y(:,ind);
        P = size(m_var,2);
        
        Rmat = m_var-(mean(m_var,2)*ones(1,P));
        
        Rmat = Rmat*Rmat'; % matrice de covariance empirique
        
        [vect_prop D] = eigs(Rmat,1);

                % Find mean of M1
        AUC = trapz(PET.time,m_var);
        
        [AUC,ind] = sort(AUC);
        quant = floor(0.1*length(ind));
        
        m_var = m_var(:,ind);
        %     SBR_M0s = double(SBR_M0s(:,2*quant:end));
        
        M_ROI = mean(m_var(:,2*quant:4*quant),2);
        close all;
end
