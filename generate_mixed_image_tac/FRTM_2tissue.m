function [tac,alpha,B] = FRTM_2tissue(vect_parameters,vect_frame_time,TAC_ref)
%% Input parameters
% vect_parameters : vecteur des parametres de la tac, dans cet ordre :
% A1 A2 A3 Lambda1 Lambda2 Lambda3 tau k1 k2 k3 k4 fv
% ex : vect_parameters=borne_sup;
%
% vect_frame_time : vecteur des points de discretisation de la tac (centre des
% frames). Le temps s'exprime en minute. Pas de dimension imposee.
% ex : t=[0.5:0.5:90];
%
% author : Yanna Cavalcanti 2018
%

%% Main

%Range des parametres si besoin
%             k1       k2       k3        k4
borne_sup = [4        0.4      0.15      0.1];
borne_inf = [1.2      0.01     0.001     0.0001];

%Variables regionales
k = vect_parameters;
R1 = k(1); %tracer delivery ratio R1 = k1/k1';
delta = sqrt((k(:,2)+k(:,3)+k(:,4)).^2-4.*k(:,2).*k(:,4));

alpha1 = (k(2)+k(3)+k(4))./2 + delta./2;
alpha2 = (k(2)+k(3)+k(4))./2 - delta./2;
B1 =  R1*(k(3) + k(4) - alpha1).*(alpha1-k(2)/R1)./delta; %B1 = b1*R1 according to armstrong
B2 =  R1*(alpha2 - k(3) - k(4)).*(alpha2-k(2)/R1)./delta; %B2 = b2*R1 according to armstrong

%Calcul de la tac
t=vect_frame_time;
exponential = (B1.*exp(-alpha1.*(t))+B2.*exp(-alpha2.*(t)))';
exponential(1) = exponential(1)+R1; 

tac = conv(exponential,TAC_ref);
tac = tac(1:length(t));

B = [R1-1,B1,B2];
alpha = [alpha1,alpha2];
