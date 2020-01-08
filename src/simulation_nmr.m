%% compare the true distribution and sampling distribution.
% Without loss of generality, we take Gaussian as example.
close all;
clear all;
clc;

%% Parameters
SampleSize = [10000,1000];
type = 'Uniform'; 

tic;
add_zero = 0;
switch type
    case 'Uniform'
        data = rand(SampleSize);
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = ones(2*N_X+1,1)/(2*N_X+1);
        % for display, add zeros for uniform
        add_zero = N_X;  
        Recover_Ex = [zeros(add_zero,1);Recover_Ex;zeros(add_zero,1)];
        Recover_Sample = [zeros(add_zero,1);Recover_Sample;zeros(add_zero,1)];
        pdf = [zeros(2*add_zero,1);pdf;zeros(2*add_zero,1)];
        N_X = N_X + add_zero;
        pdf_Ex = conv(Recover_Ex,Recover_Ex); 
    case 'Gaussian'
        data = randn(SampleSize);
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/2/pi*exp(-(i*binSize)^2/2);
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
    case 'Cauchy'
        data = tan(pi*(rand(SampleSize)-1/2));
        data = data(abs(data)<40);    % in case of too large range
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/pi/(1+(i*binSize)^2);
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
    case 'Laplacian'
        data = rand(SampleSize);
        data = sign(0.5-data).*(1/sqrt(2)).*log(2*min(data,1-data));
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/sqrt(2)*exp(-sqrt(2)*(abs(i)*binSize));
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
end
toc;
N0 = 2*N_X;

%% Display
f1 = figure;
plot([-N0:N0]*binSize,pdf,'r','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf_Ex,'b--','Linewidth',1.5);
legend([type,' From sampling'],[type,' From expression']);
set(f1,'Position',[200,300,560,420]);
grid on;

f2 = figure;
subplot(1,2,1)
plot([-N_X:N_X]*binSize,Recover_Sample,'r','Linewidth',1.5);
hold on;
plot([-N_X:N_X]*binSize,Recover_Ex,'b--','Linewidth',1.5);
legend({'$f_{\hat{X}}(t)$','$f_{X}(t)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',15);
set(f2,'Position',[760,300,1200,420]);
grid on;
title(type);
%% metric
N = length(Recover_Sample);
a = Recover_Sample;
b = Recover_Ex;
gof1.sse = sum((a-b).^2);
gof1.chisquare = nansum((a-b).^2./b);
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

subplot(1,2,2)
plot([-2*N_X:2*N_X]*binSize,conv(Recover_Sample,Recover_Sample),'r','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf,'b--','Linewidth',1.5);
legend({'$f_{\hat{X}}(t) * f_{\hat{X}}(t)$','$f_Y(t)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',15);
grid on;
title(['Autocorrelation of ',type]);
%% metric
N = length(pdf);
a = pdf;
b = pdf_Ex;
gof2.sse = sum((a-b).^2);
gof2.chisquare = nansum((a-b).^2./b);
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




