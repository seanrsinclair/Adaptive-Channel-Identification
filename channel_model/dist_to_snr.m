function SNR = dist_to_snr (Ptx,NF,d)
kb=1.38e-23; %Boltzmann cte

T=290;
BW=5e6;
Ptxlin=1e-3*10^(Ptx/10);

noiseTx=kb*T*BW;

Psignal=Ptxlin-noiseTx;

Snrin=10*log10(Psignal/noiseTx);
Snrin=50;%%tempppp

Snrout=Snrin-NF;

Snroutlin=10^(Snrout/10);

noiseTx2=(Ptxlin)/(Snroutlin+1);

noiseTx2=noiseTx2+noiseTx;

Psignal2=Ptxlin-noiseTx2;

Prx=path_loss(Ptx, 2400,d);

Prxlin=1e-3.*(10.^(Prx/10));

%decaying factor
factor=Prxlin/Ptxlin;

Psignal2=Psignal2*factor;
Pnoise=noiseTx2*factor;

% lna
Prxnoise=kb*T*BW*20;

%add noise
Pnoise=Pnoise+Prxnoise;

Snrinlna=10*log10(Psignal2./Pnoise);

%lna snr outut
Snroutlna=Snrinlna-NF;
%linear
Snroutlinlna=10.^(Snroutlna/10);

Pnoise_lna=(Psignal2+Pnoise)./(Snroutlinlna+1);
Pnoise_lna=Pnoise_lna+Pnoise;
%signal power after lna
PSignal_lna=Psignal2;

gain_lna=20;%db
gain_lna_lin=10^(gain_lna/10);

PSignal_lna=PSignal_lna*gain_lna_lin;
Pnoise_lna=Pnoise_lna*gain_lna_lin;

%mixer
Pmixernoise=0;%kb*T*BW*20;

Pnoise_mixer=Pnoise_lna+Pmixernoise;

Snrinmixer=10*log10(PSignal_lna./Pnoise_mixer);

%lna snr outut
Snroutmixer=Snrinmixer-NF;
%linear
Snroutlinmixer=10.^(Snroutmixer/10);

Pnoise_mixerout=(PSignal_lna+Pnoise_mixer)./(Snroutlinmixer+1);
Pnoise_mixerout=Pnoise_mixer+Pnoise_mixerout;
%signal power after lna
PSignal_mixer=PSignal_lna;

SNR=10*log10(PSignal_mixer./Pnoise_mixerout);

end