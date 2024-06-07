function [lleshort,llelonger]=max_lyapunov_script(x,meanperiod_samples,dt,plot_flag,save_flag)

%% init pars
sr=1/dt;
tau=round(1/4*meanperiod_samples);
maxiter=round(6*sr);
handtimevec=(1:maxiter)'.*dt;

%% div, this where the analysis forks into psr phase space, and manually reconstructed (embedded) phase space
if size(x,2)>1
    ldiv=lyarosenstein_no_psr(x,meanperiod_samples,maxiter);
else
    de=5;
    ldiv=lyarosenstein(x,de,tau,meanperiod_samples,maxiter);
end

%% fitting
indexshort=round(meanperiod_samples);
coef_short = polyfit(handtimevec(1:indexshort),ldiv(1:indexshort),1);
lleshort = coef_short(1);%*length(pseudo_time)/max(pseudo_time);

indexlonger=round(meanperiod_samples*4);
indexlongermin=round(meanperiod_samples);
coef_longer = polyfit(handtimevec(indexlongermin:indexlonger),ldiv(indexlongermin:indexlonger),1);
llelonger = coef_longer(1);%*length(pseudo_time)/max(pseudo_time);

%% visual checks
if plot_flag == 1
    if save_flag
        figure('Visible','off')
    end
    plot(handtimevec,ldiv,'-k');
    set(gcf,'DefaultTextInterpreter', 'latex')
    set(gcf,'Color','white')
    hold on
    lls=plot(handtimevec(1:indexshort),handtimevec(1:indexshort).*lleshort + coef_short(2),'r','LineWidth',2);
    lll=plot(handtimevec(indexlongermin:indexlonger),handtimevec(indexlongermin:indexlonger).*llelonger + coef_longer(2),'b','LineWidth',2);
    hold off
    xlabel('$t$ [cycles]','interpreter','latex','FontName','Arial','fontsize',16)
    ylabel('$ln(<div>)$','interpreter','latex','FontName','Arial','fontsize',16)
    set(gca,'FontName','Arial','fontsize',12)
    legend1 = legend([lls lll],['\lambda_s = ' num2str(lleshort,2)], ...
        ['\lambda_l = ' num2str(llelonger,2)], ...
        'Location','SouthEast');
    set(legend1,'LineWidth',1,'fontsize',16);
    
    if save_flag
        print('-djpeg','-r100',save_flag)
    end
end
