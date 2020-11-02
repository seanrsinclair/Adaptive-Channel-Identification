function[ Ps ] = path_loss (Ptx,fsig,d)

f = fsig; %signal frequency in MHz
Gtx = 0 ; %Tx antena gain in dBi
Grx = 10 ; %Rx antena gain in dBi

%emre: free space path loss formula
%LdB = 32.45 + 20*log10(f)+20*log10(d) -Gtx -Grx;
%emre indoor ITU-R path loss formula distance in meters
%,LdB = 20*log10(f)+30*log10(d)-28 -Gtx -Grx;
%celfree paper
%http://people.ece.cornell.edu/atang/pub/00to02/vtc2001spring.pdf
h_AP = 15;		%height of APs meter
h_u = 1.65;     %height of users meter
d0=0.01;% from the paper in km
d1=0.05;%4*h_u*h_AP/(1e-6*3e8/fsig)/1000; %in km

L = 46.3+33.9*log10(f)-13.82*log10(h_AP)-(1.1*log10(f)-0.7)*h_u+(1.56*log10(f)-0.8);

PL=zeros(length(d),1);
for i=1:length(d)
if (d(i) > d1)
    PL(i) = -L-35.*log10(d(i));
elseif (d(i) > d0)
    PL(i) = -L-15*log10(d1)-20.*log10(d(i));
else
    PL(i) = (-L-15*log10(d1)-20.*log10(d0));
end
end


Ps = Ptx +PL +Gtx + Grx;

end