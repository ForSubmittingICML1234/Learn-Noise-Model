%% Test fot the residue of signal
close all;
clear all;
clc;

%% Parameters
SNR = 10;   % 5 - 10
freq = 10;   % 1 - 20, signal changing rate
Active_ratio = 1/100;  % ratio of signal in the whole video
nonincreasing = 1;
signal_peak = 1;
N = 500*500;
T = 200;
Ns = 50/Active_ratio;
SampleSize = [N,T];


%% noise generation
noise_dif = randn(SampleSize)*sqrt(2);
[pdf_n,binSize] = getDistribution(noise_dif);
% binSize = 0.05;
% N0 = 8/binSize;
% pdf_n = zeros(2*N0+1,1);
% sigma = sqrt(2);
% for i = -N0:N0
%     pdf_n(i+N0+1) = 1/2/pi/sigma*exp(-(i*binSize)^2/2/sigma^2);
% end
% pdf_n = pdf_n/sum(pdf_n);

%% signal
sigma = 20;
interval = round(sigma/freq*5);
signal = zeros(1,2*Ns+1);
for i = -Ns:Ns
    signal(i+Ns+1) = 1/2/pi/sigma*exp(-(i)^2/2/sigma^2); % Gaussian
end
signal = signal/max(signal)*sqrt(10^(SNR/10))*sqrt(2);
signal_res = signal(1+interval:interval:end) - signal(1:interval:end-interval);
[pdf_s,binSize] = getDistribution(signal_res,binSize);
Ns2 = (length(pdf_s)-1)/2;
figure;plot([-Ns:Ns],signal,'b');hold on;plot([-Ns:interval:Ns],signal(1:interval:end),'r.','Linewidth',10);title('signal');
figure;plot([-Ns2:Ns2]*binSize,pdf_s,'Linewidth',1.5);title('pdf of signal residue');axis([-5,5,0,1]);

pdf = conv(pdf_n,pdf_s);
N0 = 1/2*(length(pdf)-1);
N1 = 1/2*(length(pdf_n)-1);
f0 = figure;
set(f0,'Position',[200,300,560*3,420]);
subplot(1,3,1);
plot([-N1:N1]*binSize,pdf_n,'b','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf,'r--','Linewidth',1.5);
legend('Without signal','With signal');
grid on;
title('pdf of Y')

%% Recover
Recover_Sample = pdf_recover(pdf);
Recover_Noise = pdf_recover(pdf_n);
N_S = (length(Recover_Sample)-1)/2;
N_N = (length(Recover_Noise)-1)/2;

subplot(1,3,2);
% set(f1,'Position',[760,300,560,420]);
plot([-N_N:N_N]*binSize,Recover_Noise,'b','Linewidth',1.5);
hold on;
plot([-N_S:N_S]*binSize,Recover_Sample,'r--','Linewidth',1.5);
legend('Without signal','With signal');
grid on;
title('Recovered pdf of X')

subplot(1,3,3);
% set(f2,'Position',[1320,300,560,420]);
plot([-2*N_S:2*N_S]*binSize,conv(Recover_Sample,Recover_Sample),'r','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf,'b--','Linewidth',1.5);
legend('Autocorrelation of X','Target Y');
grid on;
title('Comparison')

% pdf_s_rotate = [pdf_s(Ns2+1:end),zeros(1,20*(N0-Ns2)),pdf_s(1:Ns2)];
% figure;plot(real(fft(pdf_s_rotate)));

Thr = 3;
new_thr = N_S*binSize;
for i = 1:N_S
    if(sum(Recover_Sample(N_S+1+i:end))<1-normcdf(3))
        new_thr = (i-1)*binSize;
        break; 
    end
end
% end
sum(Recover_Sample(N_S+1+round(Thr/binSize):end))*100