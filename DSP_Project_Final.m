clc;
clear all;
close all;

% Import wav files:
[x,Fs]=audioread('measured_signal.wav'); 
[d,Fs]=audioread('desired_signal.wav');


% initialization
N=5000; % filter order
M=length(x); % length of x(n) and d(n)
t=1:M; 
offset=50; % NLMS offset
mu=2; % NLMS stepsize 0<mu<2
lam=0.91; % NLMS forgetting factor, also leakage
y=zeros(M,1); % initialize y(n), M*1
e=y; % initialize e(n), M*1
x1=zeros(N,1); % initialize x1(n), N*1
h=x1'; % initialize h(n), 1*N


% NLMS algorithm
for n=1:M
    x1(2:N)=x1(1:N-1); % shift temporary input signal buffer down
    x1(1)=x(n); % assign current input signal sample
    normx1=x1'*x1+offset; % update input signal vector norm
    y(n)=h*x1; % compute and assign current output signal sample
    e(n)=d(n)-y(n); % compute and assign current error signal sample
    h=h*lam+mu/normx1*e(n)*x1'; % update filter coefficient vector
end
h1=h;
MMSE=mse(e);


X=fft(x);
magX=abs(X);
D=fft(d);
magD=abs(D);
Y=fft(y);
magY=abs(Y);


% Task1 Verify cosines and their values in x(n) 
% Method 1, take fft of x(n) directly
figure;
plot(magX);
title('x(n) in frequency domain');
for n=1:20000
    s(n)=cos(2*pi*100*n/Fs)+2*cos(2*pi*200*n/Fs)+cos(2*pi*400*n/Fs)+3*cos(2*pi*500*n/Fs);
end
% Method 2, generate s(n), take xcorr of x(n) and s(n)
% s(n) is ONLY for verification
% Will have squared amplitude and exact frequencies if x(n) matches s(n)
xs=xcorr(x,s);
XS=fft(xs);
magXS=abs(XS);
figure;
plot(magXS);
title('Correlation of x(n) and s(n) in frequency domain');


% Task2 Find/estimate variance of v(n) in any way you can
% Compute variance of v via take xcorr of x and get min eigenvalue
% This part takes time if q is large!
q=5000;
r=corrmtx(x,q,'autocorrelation');
rxx=r'*r;
[V D]=eig(rxx);
alleigs=sort(diag(D));
Variance_of_v=alleigs(1)


% Task3 Find Hopt(z), its order and impulse response h(n)
% minmse=1;
% for k=1:50:10000
%     y=zeros(M,1);
%     e=y;
%     x1=zeros(k,1);
%     h=x1';
%     for n=1:M
%         x1(2:k)=x1(1:k-1); 
%         x1(1)=x(n); 
%         normx1=x1'*x1+offset; 
%         y(n)=h*x1; 
%         e(n)=d(n)-y(n);
%         h=h*lam+mu/normx1*e(n)*x1';
%     end
%     currentmse(k)=mse(e);
%     minmse=min(minmse,currentmse(k));
% end
% figure;
% plot(currentmse);
% title('MMSE with increasing filter order N');
figure;
plot(h1);
title('Impulse response h(n)');
figure;
subplot(3,1,1);
plot(t,x);
title('measured signal');
subplot(3,1,2);
plot(t,d);
title('desired signal');
subplot(3,1,3);
plot(t,y);
title('output');
figure;
subplot(3,1,1);
plot(t,magX);
title('measured signal');
subplot(3,1,2);
plot(t,magD);
title('desired signal');
subplot(3,1,3);
plot(t,magY);
title('output');


% Task4 Show how well y(n) is estimating d(n), i.e. verify and show how well amplitudes and frequencies of cosines in y(n) match those in d(n)
% Method 1, compare y(n) and d(n) directly
% time domain compare
figure;
plot(t,d,'r',t,y);
legend('desired signal','output');
axis([15000,15100,-1,1]);
title('comparison of desired signal and output');
% Method 2, determine from e(n)
figure;
plot(t,e);
title('error');
axis([1,20000,-1,1]);
SNR1=snr(x,e);
SNR2=snr(y,e);
% Method 3, compare Power Spectrum Density of y(n) and d(n)
figure;
subplot(2,1,1);
psd(d);
title('Power spectrum of desired signal');
axis([-inf,inf,-150,50]);
subplot(2,1,2);
psd(y);
title('Power spectrum of output');
axis([-inf,inf,-150,50]);


% Task5, Show the Minimum Mean Square error (MMSE) obtained for the Hopt(z)
MMSE
figure;
plot(t,d,t,y,t,e);
legend('Desired','Output','Error');
title('Desired signal, Output and Error');
xlabel('Time Index');
ylabel('Signal Value');


% play and wavwrite output signal
p=audioplayer(y,Fs); 
play(p);
wavwrite(y,'Output');