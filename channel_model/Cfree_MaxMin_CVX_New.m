clear
clc


%%  location Initiation
m = 1;		%# of APs
k = 1;			%# of Users
D = 1000;		%meter
AP_coords = D*rand(m,2);
Usr_coords = D*rand(k,2);
d = func_distance(AP_coords,Usr_coords)/1000;

%%  Parameters from Table I
N = 1;
f = 2450;		%Carrier frequency in MHz
B = 5;			%Bandwidth MHz
NF = 9;         %Noise figure dB
h_AP = 15;		%height of APs meter
h_u = 1.65;     %height of users meter
rho_d_mns = .2; %W
rho_u_mns = 0.1;
rho_p_mns = 0.1;
sgm_sh = 8;		%dB
d1 = 0.05;      %Break point 2 (km)
d0 = 0.01;      %Break point 1 (km)
tau_c = 200;    %in samples
tau_cf = 40;

%%  Noise calculation
F = 10^(NF/10);         %Noise factor
kB = 1.381e-23;         %Boltzman constant J/K
T0 = 290;               %equivalent temperature K
N_p = B*1e6*kB*T0*F;    %Hz*J/K*K so called noise power in the Journal paper
rho_d = rho_d_mns/N_p;  %so called SNR
rho_u = rho_u_mns/N_p;
rho_p = rho_p_mns/N_p;

%%  Path loss
L = 46.3+33.9*log10(f)-13.82*log10(h_AP)-(1.1*log10(f)-0.7)*h_u+(1.56*log10(f)-0.8);
PL = func_PL(d,d0,d1,L);    %dB
PL_real = 10.^(PL/10);
z = randn(m,k);
bta = PL+ sgm_sh*z;         %dB
bta_real = 10.^(bta/10);    %db to dreal domain
bta_real_real = PL_real .* 10.^(sgm_sh*z/10);   %directly calculate real domain 

%%  Channel matrix
%h = randn(m,k)+1i*randn(m,k);
h = sqrt(0.5)*(randn(m,k)+1i*randn(m,k));   %modified additive white complex Gaussian noise
g = sqrt(bta_real).*h;
