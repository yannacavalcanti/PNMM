function tac = simulation_TAC_analytique(vect_parameters,vect_frame_time)
%% Input parameters
% vect_parameters : vecteur des parametres de la tac, dans cet ordre :
% A1 A2 A3 Lambda1 Lambda2 Lambda3 tau k1 k2 k3 k4 fv
% ex : vect_parameters=borne_sup;
%
% vect_frame_time : vecteur des points de discretisation de la tac (centre des
% frames). Le temps s'exprime en minute. Pas de dimension imposee.
% ex : t=[0.5:0.5:90];
%
% author : Clovis Tauber 2012
%

%% Main 

%Range des parametres si besoin
%                A1         A2          A3         Lambda1      Lambda2        Lambda3    tau    k1       k2       k3        k4        fv                                     
    borne_sup = [20         0.8         0.4       -0.1         -0.01          -0.001      1      0.25     0.4      0.15      0.1       0.15  ];
    borne_inf = [0.1        0.01        0.01      -20          -1             -0.1        0      0.01     0.01     0.001     0.0001    0.01 ];

%Variables globales

    A1 = vect_parameters(:,1);
    A2 = vect_parameters(:,2);
    A3 = vect_parameters(:,3);
    lbd_1 = vect_parameters(:,4);
    lbd_2 = vect_parameters(:,5);
    lbd_3 = vect_parameters(:,6);
    tau = vect_parameters(:,7);
    
%Variables regionales
    k = vect_parameters(:,8:11);
    L1 = (k(:,2)+k(:,3)+k(:,4))./2 + sqrt((k(:,2)+k(:,3)+k(:,4)).^2-4.*k(:,2).*k(:,4))./2;
    L2 = (k(:,2)+k(:,3)+k(:,4))./2 - sqrt((k(:,2)+k(:,3)+k(:,4)).^2-4.*k(:,2).*k(:,4))./2;
    B1 =  k(:,1).*(k(:,3) + k(:,4) - L1)./(L2 - L1);
    B2 =  k(:,1).*(L2 - k(:,3) - k(:,4))./(L2 - L1);
    fv = vect_parameters(:,12);
    
%Calcul de la tac
    t=vect_frame_time;
    tac = zeros(size(fv));
    for i = 1:length(t);
%         if t(1,i) < tau
%             tac(1,i) = 0;
%         else
tac(:,i) = (exp(lbd_1.*(t(1,i)-tau)).*(B1./(lbd_1+L1).*(A1.*((t(1,i)-tau) - 1./(lbd_1+L1)) - A2 - A3).*(1 - fv) + B2./(lbd_1+L2).*(A1.*((t(1,i)-tau) - 1./(lbd_1+L2)) - A2 - A3).*(1 - fv) + (A1.*(t(1,i)-tau) - A2 - A3).*fv)...
            + exp(lbd_2.*(t(1,i)-tau)).*(B1.*A2./(lbd_2+L1).*(1 - fv) + B2.*A2./(lbd_2+L2).*(1 - fv) + A2.*fv)...
            + exp(lbd_3.*(t(1,i)-tau)).*(B1.*A3./(lbd_3+L1).*(1 - fv) + B2.*A3./(lbd_3+L2).*(1 - fv) + A3.*fv)...
            +exp(-L1.*(t(1,i)-tau)).*(1-fv).*(B1./(lbd_1 + L1).*(A1./(lbd_1 + L1) + A2 + A3) - A2.*B1./(lbd_2+L1) - A3.*B1./(lbd_3+L1))...
            +exp(-L2.*(t(1,i)-tau)).*(1-fv).*(B2./(lbd_1 + L2).*(A1./(lbd_1 + L2) + A2 + A3) - A2.*B2./(lbd_2+L2) - A3.*B2./(lbd_3+L2)));
%         end
    end