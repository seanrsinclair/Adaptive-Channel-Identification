function [Rx] = gen_rf_beam(s_b,fs,fc,fm,Ptx,d,N_ant,thetas)

% Written by A. Apsel (modified)
%modifications by Emre Gonultas

%d= %distance in km -emre 
%Ptx= 20; %Transmited power in dBm
% generate a qpsk signal as the input for the LNA
ts = 1/fs;

%t=0:ts:1e-3; % Time vector
% %these are not used??? (emre)
% % dfc = 1e-6*fc; % TX & RX carrier frequency mismatch in Hz
% %     dfc = 0;
% % dphi = pi*0.9; % TX & RX carrier phase mismatch in radians
% %     dphi = 0; 
% tdata = (4e-6); %bit duration 802.11(emre)
% N = 1000; %length of random binary sequence to decode
Temp=290; %Temporary temperatRATURE!!!!!!
BW=100e6;
NF=5; %noise figure (dB) 
% 
% % Baseband data generated
% % generate bpsk waveform from data. bitduration=time per bit and fc =carrier frequen
% dI =  round(rand(1,N));%generate vector of zeros and ones of length N
% dQ = round(rand(1,N)); %only use for QPSK
% 
% 
% bI=2*dI-1; % Convert unipolar to bipolar
% bQ=2*dQ-1;
% 
% sc=128;
% cp=16;
% num_ofdm_symbols = ceil(N/sc);
% 
% x_tx = zeros(sc*num_ofdm_symbols,1);
% x_tx(1:length(bI)) = bI+1j*bQ; % pad with zero symbols
% x_tx =  reshape(x_tx,sc,num_ofdm_symbols) ;
% s_o = ifft(x_tx)*sqrt(sc); % transform to time domain
% s_o(1,:)=0;
% s_o([28:36],:)=0;
% s_o = [ s_o(end-cp+1:end,:) ; s_o]; % add cyclic prefix
% bI=real(s_o(:));
% bQ=imag(s_o(:));
%%
%T=tdata; % Bit duration
t=0:(1/(fs)):2000*ts-ts; % time sequence

%%
%emre edit for pulse shaping

%imported

%%

%why not multiply with e^j2pifct ?? emre
%wI=sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier waveform I
%wQ=sqrt(2*Eb/T)*sin(2*pi*fc*t); % carrier waveform Q

%qpsk_wI=real(s_b).*wI; % modulated waveform I
%qpsk_wQ=imag(s_b).*wQ; % modulated waveform Q, note we can modify to add second data sequence later 

%this is 30dBm=1 dbW -emre
%ss = qpsk_wI+qpsk_wQ; % RF signal
%ss=(s_b.*exp(1j*2*pi*fc*t)*sqrt(2*Eb/T));
flow=2.4e9;
fhigh=2.5e9;
Rx=zeros(N_ant,length(s_b));

as=0.06;%antenna spacing [m]

dcs=zeros(N_ant,length(d));
%th=zeros(N_ant,length(d));
%th(1,:)=thetas;
%loc=[d.*sin(thetas) d.*cos(thetas)];
for ii=1:N_ant
    
    %length 
    an=(ii-1)*as;
    %distance between the next antenna and TX
    dcs(ii,:)=sqrt(d.^2+an^2-2*d*an.*cos(pi-thetas));
    %angles of each antenna
%     if ii>1
%         th(ii,:)=acos((an^2+dcs(ii,:).'.^2-d.^2)./(2*dcs(ii,:).'.*an));
%     end
    r_s=0;
    for i=1:length(fc)
        %20 dbm Ptx
        [si,siq]=mixer_beam(s_b*1e-2,50, t,fc(i),Temp,dcs(ii,i), NF);
        ss=si+1j*siq;
        %add thermal noise here at the tx side -emre
        [Prx] = path_loss(Ptx, fc(i) * 1e-6,dcs(ii,i)/1000); %Prx in dBm dcs [km]
        %received signal
        A=sqrt(1e-3*10^(Prx/10));% voltageImput 
        r_s=r_s+(A/sqrt(2))*ss; %assumes matched input
    end




    %we use a complex down converter
    [output, SNRout] = LNA(r_s,t,fm,Prx,fs,Temp,NF,BW,flow,fhigh);


    % takes the output and SNRout of LNA and passes it as the input for mixer

    %%
    [ri, rq] = mixer(output,SNRout, t,fm,Temp,fs, NF);
    Rx(ii,:)=(ri+1j*rq).';
end
%%
end