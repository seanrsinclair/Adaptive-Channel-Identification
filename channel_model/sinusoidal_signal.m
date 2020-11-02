
% generates sinusoidal signal as input for LNA

fcarrier = 920e6;
fs = 10*fcarrier;
Pin = -108; % power in dBm
t=linspace(0,100e-3,100000); % Time vector
s = cos(2*pi*fcarrier*t);
output = LNA(s,t,fcarrier,Pin);

figure(1)
plot(t, output);
xlabel('time (s)')

% takes the output of LNA and passes it as the input for mixer

[ri, rq] = mixer(output,t,fcarrier);

% filters I and Q signal

flp=1000e3;
[rbi, rbq] = low_pass_filter(ri,rq,fs,flp);

% plot baseband signal (I and Q)

figure(2)
subplot(2, 1, 1)
plot(t, ri, t, rq);
subplot(2, 1, 2)
plot(t, rbi, t, rbq);
xlabel('time (s)')


qbits=8; %number of bits
[yi,yq] = ADC(rbi,rbq,t,qbits);

figure(3)
plot(t,yi,t,yq); %axis([0 0.2e-3 -1 10])