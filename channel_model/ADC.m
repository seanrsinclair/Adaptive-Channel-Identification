
function [yi,yq]= ADC (rbi, rbq,t,qbits)

if qbits==1
    %avg=mean(signal);
    for i=1:length(rbi)
        if rbi(i)>0, yi(i)=1;
        else yi(i)=0;
        end
        if rbq(i)>0, yq(i)=1;
        else yq(i)=0;
        end
    end
else        
yi = uencode(rbi,qbits); %quantized output
yq = uencode(rbq, qbits);
end

qlevel = 1/(2^qbits);
sigma = qlevel/2;
k=2^qbits;
for j=1:k
    vdnl = normrnd(0,sigma,1,2^qbits);
    vdnl1(:,j) = vdnl';
end
vinl= sum(vdnl1');

figure(4)
plot(vinl);

%compare original data to quantized and recovered data, find BER
%berq=BER(dataQ, yi) %something hinky gong on with QPSK demod that inverts this!!!!
%beri=BER(dataI, yq)

%QBER=berq+beri %for QPSK, roughly BER


%subplot(4,1,1); 
%plot(t,bwI);axis([0 10e-6 -1 1])

%subplot(4,1,2);
%plot(t,bwQ);axis([0 10e-6 -1 1])

%subplot(4,1,3);
%plot(t,yi);axis([0 10e-6 -1 1])

%subplot(4,1,4);
% plot(t,yq);axis([0 10e-6 -1 1])

end
