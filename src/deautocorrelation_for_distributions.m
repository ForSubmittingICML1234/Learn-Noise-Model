clear all;
close all;
clc;

%% Get the histogram
% Make it as the distribution of even pdf with (2N0+1) term.
N0 = 200;
N = 2*N0+1;
pdf = zeros(2*N0+1,1);
mu = N0+1;
type = 'Linear'; 

switch type
    case 'Gaussian'
        sigma = N0/5;
        for i = 1:N
            pdf(i) = 1/2/pi/sqrt(sigma)*exp(-(i-mu)^2/2/sigma^2);
        end
    case 'Uniform'
         for i = mu-N0/2+5:mu+N0/2-5
            pdf(i) = 1;
         end     
    case 'TestA'
         for i = 1:N
            pdf(i) = 1 - 3*((i - mu)/2/N0)^2;
         end    
    case 'Rayleigh'
        sigma = 20;
         for i = 0:N-1
            pdf(i+1) = i/sigma^2*exp(-i^2/2/sigma^2);
         end    
    case 'Poisson'
        lambda = 10;
         for i = 0:N-1
            pdf(i+1) = lambda^i/gamma(i+1)*exp(-lambda);
         end   
    case 'Maxwell-Boltzmann'
        lambda = 20;
         for i = 0:N-1
            pdf(i+1) = i^2*exp(-i^2/2/lambda^2)/lambda^3;
         end 
     case 'Gamma'
        k = 2;
        theta = 20;
         for i = 0:N-1
            pdf(i+1) = 1/theta^k/gamma(k)*i^(k-1)*exp(-i/theta);
         end 
    case 'Log-normal'
        sigma = 1;
        mu0 = 1;
        ratio = 20;
         for i = 1:N
            pdf(i) = 1/i/sigma*exp(-(log(i/ratio)-mu0)^2/2/sigma^2);
         end 
    case 'Beta'
        a = 10;
        b = 100;
        ratio = N;
         for i = 1:N
            pdf(i) = 1/beta(a,b)*(i/ratio)^(a-1)*(1-i/ratio)^(b-1);
         end 
     case 'Chi-squared'
         k = 50;
         for i = 0:N-1
            pdf(i+1) = 1/2^(k/2)/gamma(k/2)*i^(k/2-1)*exp(-i/2);
         end 
    case 'Exponential'
         lambda = 1/100;
         for i = 0:N-1
            pdf(i+1) = exp(-lambda*i);
         end 
    case 'Linear'
         x0 = -50;
         for i = -N0:x0
            pdf(N0+i+1) = 1/N0/(N0+x0)*(N0+i);
         end 
         for i = x0+1:N0
            pdf(N0+i+1) = 1/N0/(N0-x0)*(N0-i); 
         end
end

pdf(isnan(pdf))=0;
pdf(isinf(pdf))=0;
pdf = pdf/sum(pdf);

%% shift center
mean_value = round([1:N]*pdf);
sum_pdf = 0;
for i = 1:N
    sum_pdf = sum_pdf + pdf(i);
    if(sum_pdf>=0.5)
        median_value = i;
        break; 
    end
end
[~,max_value] = max(pdf);
center = median_value;
pdf = [zeros(N+1-2*center,1);pdf;zeros(-N-1+2*center,1);];
if mod(length(pdf)-1,2)
   N0 = (length(pdf)-1) /2;
else
    pdf = [0;pdf;0];
    N0 = (length(pdf)-1) /2;
end
pdf_0 = pdf;

%% Recover
pdf = conv(pdf,flipud(pdf)); N0 = N0*2;
pdf = pdf/sum(pdf);

binSize = 1;
X_optimal = pdf_recover(pdf,0);
conv_X = conv(X_optimal,X_optimal);
N_X = (length(X_optimal)-1)/2;

%% Draw figures
Thr = 0.9999;
B_x = N_X*binSize;
for i = 1:N_X
   if(sum(X_optimal(N_X+1-i:N_X+1+i))>Thr)
       B_x = i*binSize;
      break; 
   end
end
B_x2 = N0*binSize;
for i = 1:N0
   if(sum(pdf(N0+1-i:N0+1+i))>Thr)
       B_x2 = i*binSize;
      break; 
   end
end

f0 = figure;
subplot(2,3,[1,4]);
plot(-N_X:N_X,X_optimal,'r','Linewidth',1.5);
hold on;
plot(-N0/2:N0/2,pdf_0,'b--','Linewidth',1.5)
legend('Recovered pdf','Target pdf'); 
axis([-B_x,B_x,0,max(X_optimal)*1.2]);
title(type);
grid on;
N = length(X_optimal);
a = X_optimal;
b = pdf_0;
gof1.sse = sum((a-b).^2);
tmp = ((a-b).^2./b);
gof1.chisquare = nansum(tmp(~isinf(tmp)));
gof1.rsquare = (N*a'*b-1).^2 / (N*sum(a.^2) -1)/(N*sum(b.^2)-1);
gof1.nrmse = sqrt(gof1.sse/N)*N;
c1 = a - 1/N;
c2 = b - 1/N;
gof1.cof = c1'*c2/sqrt(c1'*c1*c2'*c2);
ax = gca;
y_max = max(ax.YLim);
x_p = ones(1,4)*(ax.XLim(1)*0.9);
y_p = [9:-1:6]/10*y_max;
str = {{['\chi^2: ',num2str(gof1.chisquare,'%.2e')]},{['R^2: ',num2str(gof1.rsquare*100,4),'%']},{['NRMSE: ',num2str(gof1.nrmse,'%.2e')]},{['Correlation: ',num2str(gof1.cof*100,4),'%']}};
text(x_p,y_p,str);
set(f0,'Position',[200,300,560*3,420]);

subplot(2,3,[2,5]);
plot(-2*N_X:2*N_X,conv_X,'r','Linewidth',1.5);
hold on;
plot(-N0:N0,pdf,'b--','Linewidth',1.5);
legend('Autocorrelation of recovered pdf','Target Y');
title('Autocorrelation');
axis([-B_x2,B_x2,min(0,min(conv_X)*1.2),max(pdf)*1.2]);
grid on;
N = length(conv_X);
a = conv_X;
b = pdf;
gof2.sse = sum((a-b).^2);
tmp = ((a-b).^2./b);
gof2.chisquare = nansum(tmp(~isinf(tmp)));
gof2.rsquare = (N*a'*b-1).^2 / (N*sum(a.^2) -1)/(N*sum(b.^2)-1);
gof2.nrmse = sqrt(gof2.sse/N)*N;
c1 = a - 1/N;
c2 = b - 1/N;
gof2.cof = c1'*c2/sqrt(c1'*c1*c2'*c2);
ax = gca;
y_max = max(ax.YLim);
x_p = ones(1,4)*(ax.XLim(1)*0.9);
y_p = [9:-1:6]/10*y_max;
str = {{['\chi^2: ',num2str(gof2.chisquare,'%.2e')]},{['R^2: ',num2str(gof2.rsquare*100,4),'%']},{['NRMSE: ',num2str(gof2.nrmse,'%.2e')]},{['Correlation: ',num2str(gof2.cof*100,4),'%']}};
text(x_p,y_p,str);
% set(f1,'Position',[760,300,560,420]);

%% Single side
N1 = N0/2;
Single_X = zeros(N_X+1,1);
Single_Pdf = zeros(N_X+1,1);
Thr_new = nan;
for i = N_X:-1:0
    if i==N_X
        Single_X(i+1) = 2*X_optimal(N_X+1+i);
        Single_Pdf(i+1) = sum(pdf_0(1:N1+1-N_X))+sum(pdf_0(N1+1+N_X:end));
    else
        if i~=0
        Single_X(i+1) = Single_X(i+2) + 2*X_optimal(N_X+1+i);
        Single_Pdf(i+1) = Single_Pdf(i+2) + pdf_0(N1+1+i)+pdf_0(N1+1-i);
        else
            Single_X(i+1) = Single_X(i+2) + X_optimal(N_X+1);
            Single_Pdf(i+1) = Single_Pdf(i+2) + pdf_0(N1+1);
        end
    end
    
    if(isnan(Thr_new) && Single_Pdf(i+1)>2*(1-normcdf(3)))
        Thr_new = i+1;
    end
end

B_x = min(find(Single_Pdf<1e-4,1),find(Single_X<1e-4,1));
if isempty(B_x)
    B_x = N_X;
end
subplot(2,3,3);
plot(0:N_X,Single_X,'r','Linewidth',1.5);
hold on;
plot(0:N_X,Single_Pdf,'b--','Linewidth',1.5);
hold on;
plot([Thr_new,Thr_new],[0,1],'g--','Linewidth',1.5);
legend('Recovered pdf','Target pdf','0.26%');
title('Single side');
grid on;
axis([0,B_x,0,1]);
subplot(2,3,6);
ratio1 = Single_X./Single_Pdf;
ratio2 = Single_Pdf./Single_X;
ratio = max(ratio1,ratio2);
plot(0:N_X,ratio,'Linewidth',1.5);
title('Ratio');
grid on;
hold on;
plot([Thr_new,Thr_new],[0,5],'g--','Linewidth',1.5);
axis([0,B_x,0,5]);
% set(f2,'Position',[1320,300,560,420]);

% %% fft figure;
% X_optimal_rotate = [X_optimal(N0/2+1:end);X_optimal(1:N0/2)];
% pdf_rotate = [pdf_0(N0/2+1:end);pdf_0(1:N0/2)];
% figure;
% subplot(2,1,1);
% plot(real(fft(X_optimal_rotate)),'r','Linewidth',1.5);
% hold on;
% plot(real(fft(pdf_rotate)),'b--','Linewidth',1.5);
% axis([0,100,0,1]);
% grid on;
% title('Real');
% subplot(2,1,2);
% plot(imag(fft(X_optimal_rotate)),'r','Linewidth',1.5);
% hold on;
% plot(imag(fft(pdf_rotate)),'b--','Linewidth',1.5);
% axis([0,100,0,1]);
% grid on;
% title('Imag');