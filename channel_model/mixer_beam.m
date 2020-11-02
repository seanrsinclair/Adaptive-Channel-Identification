function [ri, rq] = mixer_beam(s,SNR, t,fcarrier,T,dc,NF)
%s=s+1j*imag(hilbert(s));
%input parameters

kb=1.38e-23; %Boltzmann cte


gain=8; %gain
R=50; %source impedance
IP3=10; %IIP3
P1dB=-1; %Gain compression
df=1e6;%5e6; %off set frequency for phase noise calculation
Lpn=-110; %phase noise (dBc)
if fcarrier==2.449e9
    BW=100e6;
else
BW=100e6;
noiseP=kb*T*BW;
N=length(s);
noise=wgn(1,N,noiseP,R,'linear','complex').';

hamming=0.42-0.5*cos(2*pi*[0:99]/99)+0.08*cos(4*pi*[0:99]/99);
en=norm(s);
s=conv(s,hamming,'same');

s=s/norm(s)*sqrt(2000);

s=s+noise;
end
%linearity

vip3=sqrt(2*R*10^((IP3-30)/10));
vp1db=sqrt(2*R*10^((P1dB-30)/10));
k1=gain;
k2=1;
k3=-4*k1/(3*vip3^2);
k5=-(0.1085*k1+3/4*k3*vp1db^2)*8/5/vp1db^4;
x=s;
y=k1*x+k3*(3/4)*x.^3+k5*(5/8)*x.^5;

%emre added thermal noise here(gain was added above)

%y=k1*x; %remove nonlinearity TEMPPPPPPPPPPPP1!!!!!
% Add white Gaussian noise

SNRoutmix = SNR - NF;



y = awgn(y,SNRoutmix,'measured'); %additive white Gaussian noise


%bndliit the noise at TX
if fcarrier==2.449e9
else


end
%phase noise

dphi_deg=pn_gen( fcarrier, df, Lpn);
dphi_rad=pi*dphi_deg/180;
phi_ni=normrnd(0, dphi_rad, 1,length(s));
phi_nq=normrnd(0, dphi_rad, 1,length(s));


%cri = cos(2*pi*(fcarrier)*(t-dc+l*cos(theta)/2)+phi_ni).'; %RX carrier I
%crq = sin(2*pi*(fcarrier)*(t-dc+l*cos(theta)/2)+phi_nq).'; %RX carrier Q

cri = cos(2*pi*(fcarrier)*(t-dc/3e8)+phi_ni); %RX carrier I
crq = sin(2*pi*(fcarrier)*(t-dc/3e8)+phi_nq); %RX carrier Q

rq = y.*crq.'; % ideal Demodulator Q
ri = y.*cri.';

%complex
%
%rq = imag(y.*exp(1j*2*pi*(fcarrier)*t).'); % ideal Demodulator Q
%ri = real(y.*exp(1j*2*pi*(fcarrier)*t).'); 


%subplot(2,1,1); 
%plot(t,ri);axis([0 0.05e-3 -1 1])

%subplot(2,1,2);
%plot(t,rq);axis([0 0.05e-3 -1 1])

end
