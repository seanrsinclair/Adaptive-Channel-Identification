
% Written by A. Apsel (modified)
%modifications by Emre Gonultas

% generate a qpsk signal as the input for the LNA

% Initialize parameters
close all;
clear all;
fs=100e6; % Time steps in seconds and frequency steps
ts = 1/fs;
fc = 2.495e9; %Carrier frequency in Hz
fm=2.45e9;% mixing frequency
t=0:ts:1e-3; % Time vector
%these are not used??? (emre)
% dfc = 1e-6*fc; % TX & RX carrier frequency mismatch in Hz
%     dfc = 0;
% dphi = pi*0.9; % TX & RX carrier phase mismatch in radians
%     dphi = 0; 
tdata = (4e-6); %bit duration 802.11(emre)
N = 1000; %length of random binary sequence to decode
Temp=290;
BW=100e6;

% Baseband data generated
% generate bpsk waveform from data. bitduration=time per bit and fc =carrier frequen
dI =  round(rand(1,N));%generate vector of zeros and ones of length N
dQ = round(rand(1,N)); %only use for QPSK


bI=2*dI-1; % Convert unipolar to bipolar
bQ=2*dQ-1;

sc=64;
cp=16;
num_ofdm_symbols = ceil(N/sc);

x_tx = zeros(sc*num_ofdm_symbols,1);
x_tx(1:length(bI)) = bI+1j*bQ; % pad with zero symbols
x_tx =  reshape(x_tx,sc,num_ofdm_symbols) ;
s_o = ifft(x_tx)*sqrt(sc); % transform to time domain
s_o = [ s_o(end-cp+1:end,:) ; s_o]; % add cyclic prefix
bI=real(s_o(:));
bQ=imag(s_o(:));
%%
T=tdata; % Bit duration
t=(1/(fs)):(1/(fs)):N*T/200; % time sequence

%%
%emre edit for pulse shaping
os=100;
t1=linspace(-fs,fs,os)';
s=(bI.'+1j*bQ.');
%s=[s(end:end-16) s];
on=[s; zeros(os,length(s))];
padded_array=on(:).';
ft=sinc(t1).';

s = conv( padded_array,ft);
s=s(1:2000);
%this is 30dBm=1 dbW -emre
s_b=s/norm(s)*length(s);
%%
%TX mixer
[si,siq]=mixer(s_b*1e-2,20, t,fc,Temp,fs, 20);
ss=si+1j*siq;
%add thermal noise here at the tx side -emre

Ptx= 20; %Transmited power in dBm



d=1; %distance in km -emre 
[Prx] = path_loss(Ptx, fc * 1e-6,d); %Prx in dBm

%received signal
A=sqrt(1e-3*10^(Prx/10));% voltageImput 
r_s=(A)*ss; %assumes matched input

NF=2.6; %noise figure (dB)for LNA
flow=2.4e9;
fhigh=2.5e9;
%we use a complex down converter
[output, SNRout] = LNA(r_s,t,fc,Prx,fs,Temp,NF,BW,flow,fhigh);


% takes the output and SNRout of LNA and passes it as the input for mixer
NF=10; %noise figure (dB) for mixer
[ri, rq] = mixer(output,SNRout, t,fm,Temp,fs, NF);
Rx=ri+1j*rq;
plot(fftshift(abs(fft(Rx))))