clear
clc
addpath('channel_model')
parfevalOnAll(gcp(), @warning, 0, 'off', 'all')
warning('off','all')
%rng(1
%number of channels
N=20;

%beamforming
os=100;
NN=N*os;
noise_floor=-96;%dbm

%sampling frequency
fs=100e6;
%center frequency(IF sampling)
cf=2.45e9;
f_width=5e6;
t=linspace(0,NN/fs,NN)';
%max number of active channels
sparsity=5;
%number of channels
groups=20;
kk=1;

os1=100;
t1=linspace(-fs,fs,os1)';
%data symbols
s=datasample([-1,1],NN)+1i*datasample([-1,1],NN);
on=[s; zeros(os1,length(s))];
padded_array=on(:).';
ft=sinc(t1).';
hamming=0.42-0.5*cos(2*pi*[0:99]/99)+0.08*cos(4*pi*[0:99]/99);

s1 = conv( padded_array,ft);
s=s1(1:NN).';
%data symbols-beamformed
s_b=s/norm(s)*length(s);


%%
MN=200;
M=MN/2;
SNRlist=[-5,0,5,10,15,20,25];

%distance and variance list used for average SNR calculatios
distlist=[147,95,66.5,47.3,25,13.4,1];
varlist=[135,81,46,30,30.6,18.4,9];

matload=load('dictionary_200_mutualcoherence_greedymin150.mat');
Rnew=matload.Rnew;
errorlist_aomp=zeros(3,length(distlist));
errorlist_lmp=errorlist_aomp;
errorlist_fft=errorlist_aomp;
errorlist_bomp=errorlist_aomp;
F=fft(eye(MN))/sqrt(MN);



    %channel indices in Fourier domain
    n=MN;
    inds_A=[n-(n/groups/2):n 1:n/groups/2-1 n/groups/2:n-n/groups/2-1];
    groupind=reshape(inds_A,n/groups,groups);

%list of number of wavelet coefficients  
inds_n=[100]; %100/200 = 0.5 compression
%number of monte-carlo trials
trials=500;
for ij=1:length(inds_n) % loop over subsampling factors
 
    R=Rnew(1:inds_n(ij),:);

    %Effective sensing matrix
    A=R*F';


    %%
    erroraomp=[];
    errorbomp=[];
    errorlmp=[];
    errorfft=[];
  
    y=[];

    %%
    for ii=1:length(distlist)
        error=zeros(trials,1);
        error2=zeros(trials,1);
        error3=zeros(trials,1);
        error4=zeros(trials,1);
       parfor j=1:trials
            %%
            %generate random activity
            freqs=sort(datasample(0:20-1, datasample([1:sparsity],1),'Replace',false));

            freqs_s=[freqs(freqs>=10)-10 freqs(freqs<10)+10];
            
            
            d=(distlist(ii)+rand(length(freqs),1)*varlist(ii))/1000;
                
            %RF signal
            w=gen_rf(s_b,fs,2.40e9+freqs*5e6,2.449e9,20,d).';
     
            %% a2f conversion
            x=w;
            x=x(1:MN);

            %A2F via non uniform wavelet sampling
            y(j,:)=(R*x.').';
            y2=R*w(1:MN).';

            %% Nyquist sampling result
            xf=fft(x)/sqrt(length(x));
            xis=zeros(groups,1);
            for g=1:groups
                xis(g)=norm(xf(groupind(:,g)),2);
            end
  
            [a,amin]=min(xis);
            
            %% Least Matching Pursuit
            [x_hat,aomp,lmp]=LBMP(1,4,20,A,y(j,:).');


            error(j)=length(intersect(freqs_s+1,aomp(1)));
            error2(j)=length(intersect(freqs_s+1,lmp(1)));
            error3(j)=length(intersect(freqs_s+1,amin));
            %end
        end
            target_f=zeros(trials,N);
            erroraomp(end+1)=mean(error);
            errorlmp(end+1)=mean(error2);
            errorfft(end+1)=mean(error3);
    end
    errorlist_aomp(ij,:)=erroraomp;
    errorlist_lmp(ij,:)=errorlmp;
    errorlist_fft(ij,:)=errorfft;
end
%%
for ij=1:length(inds_n)
figure
%errorxgb=[0.2452, 0.197, 0.1444, 0.1098, 0.0984, 0.0808, 0.0802, 0.0766];
semilogy(SNRlist,errorlist_aomp(ij,:))
hold on
semilogy(SNRlist,errorlist_lmp(ij,:))

semilogy(SNRlist,errorlist_fft(ij,:))
hold off
grid on
ylabel('Error rate')
xlabel('SNR')
legend('AOMP','LMP','FFT')
axis([min(SNRlist) max(SNRlist) 0.0001 0.3])
end