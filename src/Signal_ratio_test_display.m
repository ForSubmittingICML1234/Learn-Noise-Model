close all;
clear all;
clc;

%% Parameters
f0 = figure(1);
set(f0,'Position',[200,300,560,420]);
grid on;
legend_txt = cell(0);
% freqs = fliplr([1 2 5 10]);
% SNRs = [5,7.5,10,15];
Active_ratios = [1/50,1/100,1/500,1/1000];
for k = 1:numel(Active_ratios)
SNR = 10;
freq = 2;
Active_ratio = Active_ratios(k);
signal_peak = 1;
N = 500*500;
T = 200;
Ns = 50/Active_ratio;
SampleSize = [N,T];

%% noise generation
% noise_dif = randn(SampleSize)*sqrt(2);
% [pdf_n,binSize] = getDistribution(noise_dif);
binSize = 0.05;
N0 = 8/binSize;
pdf_n = zeros(2*N0+1,1);
sigma = sqrt(2);
for i = -N0:N0
    pdf_n(i+N0+1) = 1/2/pi/sigma*exp(-(i*binSize)^2/2/sigma^2);
end
pdf_n = pdf_n/sum(pdf_n);

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
pdf = conv(pdf_n,pdf_s);
Recover_Sample = pdf_recover(pdf);
N_X = (length(Recover_Sample)-1)/2;
hold on;
plot([-N_X:N_X]*binSize,Recover_Sample,'LineWidth',1.5);
legend_txt(k) = {(['Ratio: 1/',num2str(round(1/Active_ratio))])};
end

hold on;
pdf_ex = zeros(N0+1,1);
for i = -N0/2:N0/2
    pdf_ex(i+N0/2+1) = 1/2/pi*exp(-(i*binSize)^2/2);
end
pdf_ex = pdf_ex/sum(pdf_ex);
plot([-N0/2:N0/2]*binSize,pdf_ex,'LineWidth',1.5)
legend_txt(k+1) = {'True Gaussian'};
legend(legend_txt);