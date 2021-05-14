function [A0,M0,B0,alpha0,BP0] = variables_initialization(type_in,varargin)
% Define the regions for the variation study
%
%%%%%% Inputs:
% - type_in         string with the type of initialization desired
% - Y               2D image (L|N) - 'nfindr','vca','ROI nfindr','ROI vca'
% - K               number of classes - 'nfindr','vca','ROI nfindr','ROI vca'
% - M_ROI           ROI endmember - 'ROI nfindr','ROI vca'
% - m_abund         Ground truth abundances - 'ground truth','test M','test B'
% - m_eigen         Ground truth internal abundances - 'ground truth','test M','test A'
% - M_orig          Ground truth endmember (must use the transpose of it) - 'ground truth','test A','test B'
%
%%%%%% Outputs:
% - A0              Initial abundance
% - M0              Initial endmember
% - B0              Initial internal abundance
proj_simplex_array = @(y) max(bsxfun(@minus,y,max(bsxfun(@rdivide,cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);


switch type_in
    case 'nfindr'
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        A0 = NaN(K,size(PET.Y,2));
        
        %         Endmember initialization
        [M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm(PET.Y(:,PET.mask),K,'nfindr');
        % Abundance initialization (SUNSAL [2])
        A_Mask = sunsal(M0,PET.Y(:,PET.mask),'POSITIVITY','yes','ADDONE','yes');
        N = size(A_Mask,2);
        % % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
        for n = 1:N
            A_Mask(:,n) = ProjectOntoSimplex(A_Mask(:,n),1);
        end
        A0(:,PET.mask) = A_Mask;
        
        %         Internal abundance initialization
        B0 = NaN(Nv,size(PET.Y,2));
        B0(:,PET.mask) = zeros(Nv,length(PET.mask));
    case 'vca'
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        A0=NaN(K,size(PET.Y,2));
        
        %         Endmember initialization
        [M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm(PET.Y(:,PET.mask),K,'vca');
        % Abundance initialization (SUNSAL [2])
        A_Mask = sunsal(M0,PET.Y(:,PET.mask),'POSITIVITY','yes','ADDONE','yes');
        N = size(A_Mask,2);
        
        % % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
        A_Mask = proj_simplex_array(A_Mask);
        %         for n = 1:N
        %             A_Mask(:,n) = ProjectOntoSimplex(A_Mask(:,n),1);
        %         end
        A0(:,PET.mask) = A_Mask;
        
        %         Internal abundance initialization
        B0 = NaN(Nv,size(PET.Y,2));
        B0(:,PET.mask) = zeros(Nv,length(PET.mask));
    case 'ROI nfindr'
        PET = varargin{1};
        K = varargin{2};
        M_ROI = varargin{3};
        Nv = varargin{4};
        vect_prop = varargin{5};
        A0 = NaN(K,size(PET.Y,2));
        
        %         Endmember initialization
        [M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm_ROIknown(PET.Y(:,PET.mask),K,M_ROI,'nfindr');
        % Abundance initialization (SUNSAL [2])
        A_Mask = sunsal(M0,PET.Y(:,PET.mask),'POSITIVITY','yes','ADDONE','yes');
        N = size(A_Mask,2);
        % % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
        for n = 1:N
            A_Mask(:,n) = ProjectOntoSimplex(A_Mask(:,n),1);
        end
        A0(:,PET.mask) = A_Mask;
        
        %         Internal abundance initialization
        
        B0 = NaN(Nv,size(PET.Y,2));
        B0(:,PET.mask) = zeros(Nv,length(PET.mask));
    case 'ROI vca'
        PET = varargin{1};
        K = varargin{2};
        M_ROI = varargin{3};
        Nv = varargin{4};
        A0 = NaN(K,size(PET.Y,2));
        
        %         Endmember initialization
        [M0, P, U, Y_bar, endm_proj, Y_proj] = find_endm_ROIknown(PET.Y(:,PET.mask),K,M_ROI,'vca');
        % Abundance initialization (SUNSAL [2])
        A_Mask = sunsal(M0,PET.Y(:,PET.mask),'POSITIVITY','yes','ADDONE','yes');
        N = size(A_Mask,2);
        % % Abundance projection onto the unit simplex to strictly satisfy the constraints [3]
        for n = 1:N
            A_Mask(:,n) = ProjectOntoSimplex(A_Mask(:,n),1);
        end
        A0(:,PET.mask) = A_Mask;
        
        %         Internal abundance initialization
        B0 = NaN(Nv,size(PET.Y,2));
        B0(:,PET.mask) = zeros(Nv,length(PET.mask));
    case 'rand'
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        
        %         Endmember initialization
        M0 = rand(size(PET.Y,1),K);
        % Abundance initialization
        A0 = NaN(K,size(PET.Y,2));
        A_mask = rand(K,length(PET.mask));
        A0(:,PET.mask) = A_mask;
        %         Internal abundance initialization
        B0 = NaN(Nv,size(PET.Y,2));
        B_mask = rand(Nv,length(PET.mask));
        B0(:,PET.mask) = B_mask;
    case 'zero'
        
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        
        %         Endmember initialization
        M0 = zeros(size(PET.Y,1),K);
        % Abundance initialization
        A0 = NaN(K,size(PET.Y,2));
        A_mask = zeros(K,length(PET.mask));
        A0(:,PET.mask) = A_mask;
        %         Internal abundance initialization
        B0 = NaN(Nv,size(PET.Y,2));
        B_mask = zeros(Nv,length(PET.mask));
        B0(:,PET.mask) = B_mask;
    case 'ground truth'
        M_orig = varargin{1};
        m_abund = varargin{2};
        m_eigen = varargin{3};
        
        %         Endmember initialization
        M0 = M_orig';
        % Abundance initialization
        A0 = m_abund;
        %         Internal abundance initialization
        B0 = m_eigen;
    case 'test A'
        M_orig = varargin{1};
        m_eigen = varargin{2};
        
        K = size(M_orig,1);
        N = size(m_eigen,2);
        
        %         Endmember initialization
        M0 = M_orig';
        % Abundance initialization
        A0 = rand(K,N);
        %         Internal abundance initialization
        B0 = m_eigen;
    case 'test M'
        m_abund = varargin{1};
        m_eigen = varargin{2};
        L = varargin{3};
        
        K = size(M_orig,1);
        %         Endmember initialization
        M0 = rand(L,K);
        % Abundance initialization
        A0 = m_abund;
        %         Internal abundance initialization
        B0 = m_eigen;
    case 'test B'
        M_orig = varargin{1};
        m_abund = varargin{2};
        Nv = varargin{3};
        
        N = size(m_abund,2);
        %         Endmember initialization
        M0 = M_orig';
        % Abundance initialization
        A0 = m_abund;
        %         Internal abundance initialization
        B0 =zeros(Nv,N);
    case 'kmeans'
        
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        %         Endmember initialization
        M0 = zeros(size(PET.Y,1),K+2);
        % Abundance initialization
        A0 = NaN(K+2,size(PET.Y,2));
        [abund,endm] = kmeans(PET.Y(:,PET.mask)',K+2);
        for i=1:K+2
            A0(i,PET.mask) = abund==i;
            M0(:,i) = endm(i,:)';
        end
        B0 = NaN(K-1,size(PET.Y,2),Nv);
        B0(:,PET.mask,1) = eps*ones(size(B0(:,PET.mask,1)));
        B0(:,PET.mask,2) = -eps*ones(size(B0(:,PET.mask,2)));
        B0(:,PET.mask,3) = eps*ones(size(B0(:,PET.mask,3)));
        alpha0 = [0.6*ones(1,2);0.01*ones(1,2)];
        
        % ind = find(A0(1,:)==1);
        % [ B0(:,ind), Z, deltaM] = optim_B_nd(PET.Y(:,ind),A0(:,ind),M0,1,vect_prop);
    case 'segmentation'
        
        PET = varargin{1};
        K = varargin{2};
        Nv = varargin{3};
        A_gt = varargin{4};
        lbd = varargin{5};
        gamma = varargin{6};
        alpha_lim = varargin{7};
        B_lim = varargin{8};
        tracer_decay=varargin{9};
        % Abundance initialization
        A0 = NaN(size(A_gt));
        A0(:,PET.mask) = A_gt(:,PET.mask);
        
        % Endmember initialization
        M0 = zeros(size(PET.Y,1),K);
        
        % Variability initialization
        B0 = NaN(K-1,size(PET.Y,2),Nv);
        B0(:,PET.mask,1) = eps*ones(size(B0(:,PET.mask,1)));
        B0(:,PET.mask,2) = -eps*ones(size(B0(:,PET.mask,2)));
        B0(:,PET.mask,3) = eps*ones(size(B0(:,PET.mask,3)));
        alpha0 = [alpha_lim(1)*ones(1,2);alpha_lim(2)*ones(1,2)];
        for i=1:K-1
            % Find endmembers
            endm = PET.Y(:,A0(i,:)==1);
            AUC = trapz(PET.time,endm);
            [AUC,ind] = sort(AUC);
            quant = floor(0.1*length(ind));
            endm = endm(:,ind);
            M0(:,i) = mean(endm(:,quant:2*quant),2);
        end
        M0(:,K) = mean(PET.Y(:,A0(K,:)==1),2);
        [B0,alpha0,BP0,Q0,elapsedTime,f0,i0] = basis_fcts_init(PET,A0,M0,B0,alpha0,3000,lbd,gamma,alpha_lim,B_lim,1e-04,tracer_decay);

% Remember pe2i
% 1e-5
% err0 =
%     0.4255    0.1475    2.0375    0.2780
% 1e-4
% err0 =
%          0    0.1475    1.2523    0.1968
% 1e-3
% err0 =
%     0.2789    0.1475    0.8411    0.0631
% 1e-2
% err0 =
%     0.2789    0.1475    0.8145    0.0970
end