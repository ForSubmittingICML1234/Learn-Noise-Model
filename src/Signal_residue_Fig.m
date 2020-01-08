%% Test fot the residue of signal
close all;
clear all;
clc;

%% Parameters
SNRs = [5,5,5,5,5,7.5,10,10];   % 5 - 10
freqs = [20,10,5,5,5,5,5,3];   % 1 - 20, signal changing rate
Active_ratios = [200,200,200,100,50,50,50,50];  % ratio of signal in the whole video
SampleSize = 10e8;
FalsePostiveRates = [5e-5,7.5e-5,1e-4,2.5e-4,5e-4,7.5e-4,1e-3,2.5e-3,5e-3,7.5e-3,1e-2,2.5e-2,5e-2];

nonincreasing = 1;
signal_peak = 1;



%% noise generation
% [pdf_n,binSize] = getDistribution2(SampleSize);
% save('Gaussian_Noise_1e9','pdf_n','binSize');
load('Gaussian_Noise_1e9.mat');
% clear noise_dif;
% results = cell(0);
%% signal
figure;
loglog(FalsePostiveRates,FalsePostiveRates,'Linewidth',1.5,'DisplayName','Standard Gaussian');
axis([1e-4,1e-1,1e-6,1e-1]);

xlabel({'$\int_t^{\infty}f_{\hat{N}}(x)dx$'},'Interpreter','latex','FontName', 'Calibri','FontSize',12);
ylabel({'$\int_t^{\infty}f_{{N}}(x)dx$'},'Interpreter','latex','FontName', 'Calibri','FontSize',12);
% NN = (length(pdf_n)-1)/4;
% pdf = zeros(NN*2+1,1);
% for i = 1:numel(pdf)
%     pdf(i) = exp(-((i-1-NN)*binSize)^2/2);
% end
% pdf = pdf / sum(pdf);
% % pdf_n = conv(pdf,pdf);
% plot([-NN:NN]*binSize,pdf,'Linewidth',1.5,'DisplayName','Standard Gaussian');
% axis([2.5,5.5,0,1e-3]);
hold on;
for kkk = 1:numel(SNRs)
    SNR = SNRs(kkk);
    freq = freqs(kkk);
    Active_ratio = 1/Active_ratios(kkk);
    Ns = 60/Active_ratio;
    sigma = 20;
    interval = round(sigma/freq*6);
    signal = zeros(1,2*Ns+1);
    for i = -Ns:Ns
        signal(i+Ns+1) = 1/2/pi/sigma*exp(-(i)^2/2/sigma^2); % Gaussian
    end
    signal = signal/max(signal)*sqrt(10^(SNR/10))*sqrt(2);
    signal_res = signal(1+interval:interval:end) - signal(1:interval:end-interval);
    [pdf_s,binSize] = getDistribution(signal_res,binSize);
    Ns2 = (length(pdf_s)-1)/2;
%     figure;plot([-Ns:Ns],signal,'b');hold on;plot([-Ns:interval:Ns],signal(1:interval:end),'r.','Linewidth',10);title('signal');
    % figure;plot([-Ns2:Ns2]*binSize,pdf_s,'Linewidth',1.5);title('pdf of signal residue');axis([-5,5,0,1]);

    pdf = conv(pdf_n,pdf_s);
    %% Recover
    Recover_Sample = pdf_recover(pdf);
    N_S = (length(Recover_Sample)-1)/2;
    Real_FPRs = zeros(numel(FalsePostiveRates),1);
    for k = 1:numel(FalsePostiveRates)
        FalsePostiveRate = FalsePostiveRates(k);
        for i = N_S:-1:1
            if(sum(Recover_Sample(N_S+1+i:end))>FalsePostiveRate)
                thr = i*binSize;
                Real_FPRs(k) = 1-normcdf(thr);
                break; 
            end
        end
    end
    
    
    name = ['Freq ', num2str(freq),', SNR ',num2str(SNR),'dB, Ratio 1/',num2str(Active_ratios(kkk))];
    loglog(FalsePostiveRates,Real_FPRs,'Linewidth',1.5,'DisplayName',name);
%     plot([-N_S:N_S]*binSize,Recover_Sample,'Linewidth',1.5,'DisplayName',name);
%     hold on;
end
grid on;
legend;