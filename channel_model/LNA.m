
function [ output, SNRout ] = LNA( x, t, fcarrier, Pin,fs,T,NF,BW,flow,fhigh)

%receives an input signal s with power Pin in dBm, and carrier frequency
%fcarrier.

%input parameters

kb=1.38e-23; %Boltzmann cte
gain=20; % voltage gain (dB)
R=50; %source impedance
IP3=4.7; %IIP3 (dBm)
P1dB=-1; %Gain compression (dBm)
%BW=1e9;




% nonlinearity

vip3=sqrt(2*R*10^((IP3-30)/10));
vp1db=sqrt(2*R*10^((P1dB-30)/10));
gain =10^(gain/20);
k1=gain;
k2=0;
k3=-4*k1/(3*vip3^2);
k5=-(0.1085*k1+3/4*k3*vp1db^2)*8/5/vp1db^4;
y=k1*x+k2*x.^2+k3*(3/4)*x.^3+k5*(5/8)*x.^5;
%y=k1*x; %remove nonlinearity TEMPPPPPPPPPPPP1!!!!!

%noise

noiseP=kb*T*BW*gain^2*(10^(NF/10)-1);
N=length(t);
noise=wgn(1,N,noiseP,R,'linear','complex').';

%bln= band_pass_filter( noise, t, flow, fhigh, fcarrier*10 );%band limits the noise
%noise=(sqrt(noiseP)/std(bln))*bln;%rescale band limited noise to reflect NF correctly
%output= band_pass_filter( y, t, 1, 5e6, fs );%band limits the signal and noise
output=y+noise;

% input and output SNR 

SNRin=1e-3*10^(Pin/10)/(BW*kb*T*10^(NF/10));
SNRout= SNRin/(10^(NF/10));
SNRin=10*log10(SNRin);
SNRout=10*log10(SNRout);


%subplot(4,1,1);
%plot(t,s,t,y);

%subplot(4,1,2);
%plot(t,y);

%subplot(4,1,3);
%plot(t,noise);

%subplot(4,1,4);
%plot(t,output);
%xlabel('time (s)')
end