%% Creation of binary wavelets
clear

N=200;  %Length of the input signal

%generate mother wave
Wave=zeros(1,N);
Wave(1:2:end)=1;
Wave=2*Wave-1;
ll=0;
A=[];
%% NUWS
Mother=[];
divisors=[1; unique(cumprod(perms(factor(N)),2))];
f=divisors(1:end-1);
t=linspace(0,1,N);
%generate mother waves with different frquencies
for l=f
    Mother=[Mother ;2*(sin(2*pi*l*t)>=0)-1];
end
A=[]
scales=[f];
for i=1:length(scales)
    A=[A;kron(eye(scales(i)),Mother(i:end,1:N/scales(i)))];
end
A=unique(A,'rows','stable');
A=(A(sum(A,2)==0,:))
imagesc(A)