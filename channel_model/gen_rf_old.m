function [Rx] = gen_rf_old(fs,fc,fm,d)

% Written by A. Apsel (modified)
%modifications by Emre Gonultas

%d= %distance in km -emre 
Ptx= 20; %Transmited power in dBm
% generate a qpsk signal as the input for the LNA
ts = 1/fs;

t=0:ts:1e-3; % Time vector
%these are not used??? (emre)
% dfc = 1e-6*fc; % TX & RX carrier frequency mismatch in Hz
%     dfc = 0;
% dphi = pi*0.9; % TX & RX carrier phase mismatch in radians
%     dphi = 0; 
tdata = (4e-6); %bit duration 802.11(emre)
N = 5; %length of random binary sequence to decode
Temp=290;
BW=100e6;

% Baseband data generated
% generate bpsk waveform from data. bitduration=time per bit and fc =carrier frequen
dI =  round(rand(1,N));%generate vector of zeros and ones of length N
dQ = round(rand(1,N)); %only use for QPSK


bI=2*dI-1; % Convert unipolar to bipolar
bQ=2*dQ-1;

T=tdata; % Bit duration
Eb=T/2; % This will result in unit amplitude waveforms

t=(1/(fs)):(1/(fs)):N*T; % time sequence
Ns=length(t); % Number of samples
Nsb=Ns/N; % Number of samples per bit
ddI=repmat(dI',1,Nsb); % replicate each bit Nsb times
ddQ=repmat(dQ',1,Nsb);
bbI=repmat(bI',1,Nsb); dwI=ddI'; % Transpose the rows and columns
bbQ=repmat(bQ',1,Nsb); dwQ=ddQ';
dwI=dwI(:)'; 
dwQ=dwQ(:)';
% Convert dw to a column vector (colum by column) and convert to a row vector
bwI=bbI';
bwQ=bbQ';

bwI=bwI(:)'; % Data sequence samples
bwQ=bwQ(:)';
%%
%emre edit for pulse shaping
os=100;
t1=linspace(-fs,fs,os)';
s=bwI+1j*bwQ;
on=[s; zeros(os,length(s))];
padded_array=on(:).';
ft=sinc(t1).';

s = conv( padded_array,ft);
s=s(1:Ns);
s_b=s/norm(s)*length(s);
%%

%why not multiply with e^j2pifct ?? emre
%wI=sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier waveform I
%wQ=sqrt(2*Eb/T)*sin(2*pi*fc*t); % carrier waveform Q

%qpsk_wI=real(s_b).*wI; % modulated waveform I
%qpsk_wQ=imag(s_b).*wQ; % modulated waveform Q, note we can modify to add second data sequence later 

%this is 30dBm=1 dbW -emre
%ss = qpsk_wI+qpsk_wQ; % RF signal
%ss=(s_b.*exp(1j*2*pi*fc*t)*sqrt(2*Eb/T));

%20 dbm Ptx
[si,siq]=mixer(s_b*1e-2,30, t,fc,Temp,fs, 10);
ss=si+1j*siq;
%add thermal noise here at the tx side -emre



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

end