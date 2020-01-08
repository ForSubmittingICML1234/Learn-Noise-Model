%% For figure 2
close all;
clear all;
clc;

%% Parameters
SampleSizes = [502*502*103,502*502*1614,502*502*180];
figNum = 3;

%% Display
f = figure;
set(f,'Position',[0,0,1400,1080]);
for k = 1:figNum
    filename = ['pdf',num2str(k),'.mat'];
    load(filename);
    weightType = 'Average';
    Recover_Sample = pdf_recover(pdf,SampleSizes(k),weightType,0);
    N_X = (length(Recover_Sample)-1)/2;
    pdf_Ex = conv(Recover_Sample,Recover_Sample);
    N0 = (length(pdf)-1)/2;
    

    subplot(2,figNum,k,'Position',[-0.28+k*0.32 0.56 0.27 0.38])
    plot([-N_X:N_X]*binSize,Recover_Sample,'r','Linewidth',1.5);
    legend({'$f_{\hat{N}}(x)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',12);
    grid on;
    title(['Dataset ',num2str(k)]);

    subplot(2,figNum,figNum+k,'Position',[-0.28+k*0.32 0.07 0.27 0.4])
    plot([-2*N_X:2*N_X]*binSize,pdf_Ex,'r','Linewidth',1.5);
    hold on;
    plot([-N0:N0]*binSize,pdf,'b--','Linewidth',1.5);
    legend({'$f_{\hat{N}}(x) * f_{\hat{N}}(x)$','$f_D(x)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',12);
    grid on;
    title(['Autocorrelation']);
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

end


