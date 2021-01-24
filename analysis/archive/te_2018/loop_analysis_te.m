function S=loop_analysis_te(flbase,fl,OutputFolder)

plotting=0;
fn=dir([flbase '/' fl 'movt*']);
S =nan(length(fn),16);
pp=str2double(fl(regexpi(fl,'pp')+2:regexpi(fl,'pp')+3));
condition=str2double(fl(regexpi(fl,'group')+5:regexpi(fl,'group')+6));

trial_counter=zeros(1,3);
for k=1:length(fn)
    fprintf('%10s,%70s\n',fl,fn(k).name)
    trial=k;
    pre=~isempty(regexpi(fn(k).name,'pre'));
    post=~isempty(regexpi(fn(k).name,'post'));
    stim_type=str2double(fn(k).name(end));
    if pre==1
        type=str2double(fn(k).name(regexpi(fn(k).name,'_pre_p')+6));
    elseif post==1
        type=str2double(fn(k).name(regexpi(fn(k).name,'_post_p')+7));
    else
        type=0;
    end
    epsilon=0;
    if all([~pre,~post,condition==3])
        epsilon=str2double(fn(k).name(regexpi(fn(k).name,'eps')+(4:6)));
    end
    x=dlmread([flbase '/' fl fn(k).name],',',1,0);
    x=x(x(:,1)>10,:);
    
    if any([pre post])
        if stim_type==3
            stim_type = 2;
            m_stim=3;
            m_pp=3;
        else
            stim_type = 1;
            m_stim=2;
            m_pp=2;
        end
    else
        stim_type = 3;
        m_stim=3;
        m_pp=3;
    end
    
    tau_stim=quarter_period(x,11);
    tau_pp=quarter_period(x,12);
    
    S(k,1:11)=[pp,trial,condition,pre,post,type,epsilon,tau_stim,m_stim,tau_pp,m_pp];
    S(k,12)=rms(x(:,11)-x(:,12));
    %r=sync(x(:,11),x(:,12),tau,m,nn,theiler);
    S(k,13:16) = 0;

    trial_counter(stim_type)=trial_counter(stim_type)+1;
    DATA{stim_type}.data.trial{trial_counter(stim_type)}=x(:,[11 12])';
    DATA{stim_type}.data.time{trial_counter(stim_type)}=x(:,1)'-x(1,1);
    DATA{stim_type}.data.fsample(trial_counter(stim_type))=1/mean(diff(x(:,1)));
    DATA{stim_type}.data.label={'Tutor','PP'};
    DATA{stim_type}.data.tau{trial_counter(stim_type)}=[tau_stim;tau_pp];
    DATA{stim_type}.data.dim{trial_counter(stim_type)}=[m_stim;m_pp];
end

for stim_type=1:3
    fileoutsave=fullfile(OutputFolder,['RawDataForTrentool_' fl(1:end-1) '_stim_type_' num2str(stim_type)]);
    data=DATA{stim_type}.data;
    save(fileoutsave,'data')
end

fprintf('\n');

if 0
    labels2=[{'S(CPG|User)'},'S(User|CPG)','H(CPG|User)','H(User|CPG)','N(CPG|User)','N(User|CPG)','CC'];
    labels=[{'pp'},'trial','condition','pre','post','type','\epsilon','\langleCrossCorr\rangle','RMS_{X-Y}','C_{max}/RMS_{X-Y}'];
    %     if condition==1
    %         try scores_training=dlmread([flbase fl 'scores_Lorenz.csv'],',',1,0);catch; scores_training=dlmread([flbase fl 'scores_Sine.csv'],',',1,0);end
    %     else
    %         scores_training=dlmread([flbase fl 'scores_Chua.csv'],',',1,0);
    %     end
    %     if size(scores_training,2)==14;
    %         scores_training=scores_training(:,2:14);
    %     end
    
    %try xs=dlmread([flbase fl 'scores_sample_entropy.csv'],',',1,0);catch;xs=dlmread([flbase fl 'scores_sample_entropy.csv'],';',1,0);end;
    %if mean(xs(:,1))==0;xs(:,1)=[];end
    %xs=[xs(xs(:,4)==1,:);xs(sum(xs(:,4:5),2)==0,:);xs(xs(:,5)==1,:)];
    %%S(k,8)=ccmax;
    %S(:,8)=xs(:,8);
    %%S(k,10)=ccmax/rms;
    %S(:,10)=xs(:,16);
    
    labels=[labels labels2];
    fid=fopen([flbase 'scores_nl_inter.csv'],'w');
    for kk=1:length(labels)
        fprintf(fid,'%20s;',labels{kk});
    end
    fprintf(fid,'\n');
    for kk=1:size(S,1)
        fprintf(fid,'%20.4f;',S(kk,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end

if plotting==1
    figure(1)
    
    subplot(1,2,1);
    scatter(S(:,9),S(:,end-4));
    title('H(CPG|User) ~ RMSE')
    xlabel('RMSE');ylabel('H(CPG|User)')
    
    subplot(1,2,2);
    scatter(S(:,9),S(:,end-3));
    title('H(User|CPG) ~ RMSE')
    xlabel('RMSE');ylabel('H(User|CPG)')
    
    pause
end
end

function tau=quarter_period(x,col)
dt=mean(diff(x(:,1)));
if isempty(col)
    if mean(diff(x(:,8)))==0
        xx=x(:,9)-mean(x(:,9));
    else
        xx=x(:,8)-mean(x(:,8));
    end
else
    xx=x(:,col)-mean(x(:,col));
end
%xx=abs(xx);
xx=detrend(xx);
xx=xx./max(xx);
[~, pklocs_cpg]=findpeaks(xx,'MinPeakHeight',.2,'MinPeakDistance',.5*1/dt);

% plot(xx);hold on;plot(pklocs_cpg,xx(pklocs_cpg),'or');hold off;pause

tau=round(mean(diff(pklocs_cpg))/4);
end