function S=loop_analysis_lyap_divergence(flbase,fl,plotting,save_flag)

filt_freq=6;
trial_transient = 10;

fn=dir(fullfile(flbase,fl,'trial_log*.txt'));
S=nan(length(fn),8);

parfor k=1:length(fn)
%for k=1:length(fn)
    pp = str2double(fl(regexpi(fl,'pp')+(2:3)));
    condition = str2double(fn(k).name(regexpi(fn(k).name,'task')+(5:6)));
    trial = k;
    x = dlmread(fullfile(flbase,fl,fn(k).name),',',1,0);
    x = x(x(:,1) > trial_transient,:);
    dt = mean(diff(x(:,1)));
    
    [l,c]=xcov(x(:,11),200,'coef');
    [~,pksloc]=findpeaks(l(c>20)-mean(l(c>20)),'MinPeakHeight',.1,'MinPeakDistance',50);
    if plotting==2
        figure(2)
        plot(c,l)
        hold on
        plot(c(pksloc + (numel(c)-1)/2+1+10),...
            l(pksloc + (numel(c)-1)/2+1+10),'v')
        hold off
    end
    [~,pksloc]=findpeaks(l(c>20)-mean(l(c>20)),'MinPeakHeight',.1,'MinPeakDistance',50);
    pksloc = pksloc(1);
    
    meanperiod = c(pksloc + (numel(c)-1)/2+1+10)*dt*1e3;
    
    [bf,af] = butter(4,filt_freq/(1/dt/2));
    % x(:,7:9) = filtfilt(bf,af,x(:,7:9));
    x(:,2:4) = filtfilt(bf,af,x(:,2:4));
    x(:,6) = filtfilt(bf,af,x(:,6));
    
    fprintf('%s, %s\n','tutor',fn(k).name)
    if save_flag == 1
        save_name = [num2str(pp,'%02.0f') '_' fn(k).name(11:end-5) '_tutor_' datestr(now,'yyyymmddHHMMSS') '.jpg'];
    else
        save_name = false;
    end
    if condition==10
        s = -cos(x(:,7));
        v = diff(s)./dt;
        v = filtfilt(bf,af,v);
        a = diff(v)./dt;
        a = filtfilt(bf,af,a);
        [lyaptutor,lyaptutorl] = max_lyapunov_script([s(3:end),v(2:end),a(1:end)],round(meanperiod/dt/1e3),dt,plotting,save_name);
    else
        [lyaptutor,lyaptutorl] = max_lyapunov_script(x(:,7:9),round(meanperiod/dt/1e3),dt,plotting,save_name);
    end
    fprintf('\n%s, %s\n','learner',fn(k).name)
    if save_flag == 1
        save_name = [num2str(pp,'%02.0f') '_' fn(k).name(11:end-5) '_learner_' datestr(now,'yyyymmddHHMMSS') '.jpg'];       
    else
        save_name = false;
    end
    [lyaplearners,lyaplearnerl]=max_lyapunov_script(x(:,6),round(meanperiod/dt/1e3),dt,plotting,save_name);
    S(k,:)=[pp trial condition meanperiod lyaptutor lyaptutorl lyaplearners lyaplearnerl];
end
fprintf('\n');

if plotting==3
    figure(3)
    
    subplot(2,2,1)
    scatter(S(:,2),S(:,end-1))
    title('$\lambda_{short,Trainee}$','interpreter','latex')
    xlabel('Trial')
    ylabel('$\langle log(divergence) \rangle$')
    
    subplot(2,2,2)
    scatter(S(:,2),S(:,end-0))
    title('$\lambda_{long,Trainee}$')
    xlabel('Trial')
    ylabel('$\langle log(divergence) \rangle$')
    
    subplot(2,2,3)
    scatter(S(:,2),S(:,end-3))
    title('$\lambda_{short~Tutor}$')
    xlabel('Trial')
    ylabel('$\langle log(divergence) \rangle$')
    
    subplot(2,2,4)
    scatter(S(:,2),S(:,end-2))
    title('$\lambda_{long~Trainee}$')
    xlabel('Trial')
    ylabel('$\langle log(divergence) \rangle$')
end
