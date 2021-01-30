%% set software paths
addpath('~/Desktop/TRENTOOL3')
addpath('~/Desktop/TRENTOOL3/fieldtrip');
addpath('~/Desktop/umcg_hidden_dynamics/analysis/te_2021/')
ft_defaults;

%% define data paths
base_folder = '/home/lilabuntu/Desktop/umcg_data';
OutputDataPath = fullfile(base_folder,'results/');
InputDataPath  = fullfile(base_folder,'training_processed/training');

%load(fullfile(InputDataPath,'lorenz_1-2_45ms.mat'));
%DATA = convert_csv_to_ft_and_te(InputDataPath);
cd(fullfile(base_folder))

TGA_results = cell(1,numel(DATA));
parfor pp = 1:numel(TGA_results)
%for pp = 1%:numel(TGA_results)
    %% data
    data = DATA{pp};
    [cfgTEP,cfgTESS] = define_cfgs(data);
    cfgTEP.verbosity = 'none';
    cfgTESS.optdimusage = 'maxdim';
    cfgTESS.numpermutation = 'findDelay';
    
    % results file name
    cfgTESS.fileidout  = strcat(OutputDataPath,'out_');
    
    %% TE calculation - scan over specified values for u
    TGA_results{pp} = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);
    ntrials = size(TGA_results{pp}.TEmat,2);
    ncomps = size(TGA_results{pp}.TEmat,1);
    maxnsurrs = ntrials^2-ntrials;
    TGA_results{pp}.TEsurr = nan(ncomps,ntrials,ntrials);
    TGA_results{pp}.TEsurr_prop = nan(size(TGA_results{pp}.TEmat));
    
    %% Exhaustive trial-matching surrogates
    cfgTESS.trial_shuffle_surr_n_per_trial = 10;
    nsurrs = ntrials*cfgTESS.trial_shuffle_surr_n_per_trial;
    counter = 0;
    for tr1 = 1:ntrials
        trials2 = 1:ntrials;
        trials2(trials2==tr1) = [];
        trials2 = trials2(randperm(numel(trials2),cfgTESS.trial_shuffle_surr_n_per_trial)');
        for tr2ind = 1:numel(trials2)
            tr2 = trials2(tr2ind);
            if tr1~=tr2
                counter = counter + 1;
                
                datatr = data;
                datatr.trial = data.trial(tr1);
                %datatr.trial{1}(1,:) = data.trial{tr1}(1,:);
                datatr.trial{1}(2,:) = data.trial{tr2}(2,:);
                datatr.time = datatr.time(tr1);
                cfgTEPtr = cfgTEP;
                cfgTEPtr.predicttimemin_u    = 10;      % minimum u to be scanned
                cfgTEPtr.predicttimemax_u    = 1010;	  % maximum u to be scanned
                cfgTEPtr.predicttimestepsize = 50; 	  % time steps between u's to be scanned
                cfgTEPtr.ragdim = max(TGA_results{pp}.cfg.dim);
                cfgTEPtr.minnrtrials = 0;   % minimum acceptable number of trials
                TGAtemp = InteractionDelayReconstruction_calculate(cfgTEPtr,cfgTESS,datatr);
                TGA_results{pp}.TEsurr(:,tr1,tr2) = TGAtemp.TEmat;
                if mod(counter,10)==0;fprintf('%6.2f%% surrs done.\n',counter/nsurrs*1e2);end
            end
        end
    end
    surrs = reshape(TGA_results{pp}.TEsurr,2,[],1);
    for tr = 1:ntrials
        TGA_results{pp}.TEsurr_prop(:,tr) = (sum(TGA_results{pp}.TEmat(:,tr)'>surrs')./nsurrs)';
    end
end
%save([OutputDataPath 'out_TGA_results.mat'],'TGA_results','DATA','OutputDataPath','InputDataPath');
%save([OutputDataPath 'out_TGA_results.mat'],'TGA_results','DATA','OutputDataPath','InputDataPath','TEtable');


return

%% Collect in long table for stats.
TEtable = [];
for pp=1:numel(TGA_results)
    TEtable = vertcat(TEtable,[DATA{pp}.conditions,...
        (1:size(DATA{pp}.conditions,1))',...
        TGA_results{pp}.TEmat', ...
        DATA{pp}.dim', ...
        DATA{pp}.tau'./DATA{pp}.fsample, ...
        DATA{pp}.scores, ...
        TGA_results{pp}.cfg.dim'.*ones(size(TGA_results{pp}.TEmat,2),1), ...
        (squeeze(TGA_results{pp}.ACT.actvalue(1,:,:))'.*TGA_results{pp}.cfg.tau)./DATA{pp}.fsample, ...
        TGA_results{pp}.TEpermvalues(:,1)'.*ones(size(TGA_results{pp}.TEmat,2),1), ...
        TGA_results{pp}.TEpermvalues(:,6)'.*ones(size(TGA_results{pp}.TEmat,2),1), ...
        TGA_results{pp}.TEsurr_prop']);
end
column_labels = {'Date','Time','Task','Auditory','Visual','epsilon',...
    'trial',...
    'TEstim','TEpp',...
    'dimstim','dimpp','taustim','taupp',...
    'C','rmse','tau','cbyrmse',...
    'dimhatstim','dimhatpp','tauhatstim','tauhatpp',...
    'permtest12','permtest21',...
    'interactdelay12','interactdelay21','surrprop12','surrprop21'};
TEtable = array2table(TEtable,'VariableNames',column_labels);
TEtable.Task(TEtable.Task==10)=1;
TEtable.Task(TEtable.Task==30)=2;
TEtable.Task(TEtable.Task==25)=3;
%{
boxplot(TEtable.interactdelay12,TEtable.Task)
boxplot(TEtable.interactdelay21,TEtable.Task)
boxplot(TEtable.dimhatpp,TEtable.Task)
boxplot(TEtable.dimhatstim,TEtable.Task)
boxplot(TEtable.taupp,TEtable.Task)
boxplot(TEtable.tauhatpp,TEtable.Task)
boxplot(TEtable.TEstim,TEtable.Task)
boxplot(TEtable.TEpp,TEtable.Task)
boxplot(TEtable.surrprop12,TEtable.Task)
boxplot(TEtable.surrprop21,TEtable.Task)
boxplot(TEtable.permtest12,TEtable.Task)
boxplot(TEtable.permtest21,TEtable.Task)
%}


%% Visually check the TE and surrogate distibutions and stats.
for pp=1:numel(TGA_results)
    if isfield(TGA_results{pp},'TEmat_sur')
        for d=1:size(TGA_results{pp}.TEmat_sur,1)
            [c0,n0]=hist(TGA_results{pp}.TEmat_sur(d,:));
            [c,n]=hist(TGA_results{pp}.TEmat(d,:));
            subplot(3,1,1)
            plot(n0,c0,'--',n,c,'-')
            hold on
        end
        hold off
        legend('stim0->user0','stim->user','user0->stim0','user->stim')
    end
    subplot(3,1,2)
    plot(TGA_results{pp}.TEmat')
    subplot(3,1,3)
    plot(TGA_results{pp}.TEsurr_prop')
    fprintf('%8.0f%8.0f%8.0f%8.0f%8.0f%8.3f\n',mean(DATA{pp}.conditions))
    fprintf('\n')
    disp(TGA_results{pp}.TEpermvalues)
    disp(TGA_results{pp}.cfg.dim')
    pause
end
% scatter(TEtable.trial(TEtable.Task==3),TEtable.surrprop21(TEtable.Task==3));


%% Check trends. Average TE across pps per training trial.
TE = zeros(40,2,3);
TEcount= zeros(40,1,3);
for f=1:numel(DATA)
    switch DATA{f}.conditions(1,3)
        case 10
            task = 1;
        case 30
            task = 2;
        case 25
            task = 3;
    end
    TEcount(1:numel(TGA_results{f}.TEmat(1,:)'),1,task) = TEcount(1:numel(TGA_results{f}.TEmat(1,:)'),1,task)+1;
    TE(1:numel(TGA_results{f}.TEmat(1,:)'),1,task) = TE(1:numel(TGA_results{f}.TEmat(1,:)'),1,task)+TGA_results{f}.TEmat(1,:)';
    TE(1:numel(TGA_results{f}.TEmat(1,:)'),2,task) = TE(1:numel(TGA_results{f}.TEmat(1,:)'),2,task)+TGA_results{f}.TEmat(2,:)';
    %TE(1:numel(TGA_results{f}.TEmat(1,:)'),1,task) = TE(1:numel(TGA_results{f}.TEmat(1,:)'),1,task)+TGA_results{f}.TEsurr_prop(1,:)';
    %TE(1:numel(TGA_results{f}.TEmat(1,:)'),2,task) = TE(1:numel(TGA_results{f}.TEmat(1,:)'),2,task)+TGA_results{f}.TEsurr_prop(2,:)';
end
for task=1:3
    TE(:,1,task) = TE(:,1,task)./TEcount(:,1,task);
    TE(:,2,task) = TE(:,2,task)./TEcount(:,1,task);
    TE(:,3,task) = task;
end
TE = reshape(permute(TE,[1 3 2]),[],3);
for task=1:3
    subplot(1,3,task)
    plot(TE(TE(:,3)==task,1),'-ob')
    hold on
    plot(TE(TE(:,3)==task,2),'-or')
    plot(abs(diff(TE(TE(:,3)==task,1:2),1,2)),'-dm')
    hold off
    legend('stim->user','user->stim','abs(stim-user)')
end

