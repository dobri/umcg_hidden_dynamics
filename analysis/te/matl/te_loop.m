clear

%% set software paths
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/')
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/libs/TRENTOOL3')
addpath('~/logos/umcg_hidden_dynamics/analysis/te/matl/libs/fieldtrip');
ft_defaults;

%% define data paths
base_folder = '~/biomech/projects/side_projects/umcg/data/';
OutputDataPath = fullfile(base_folder,'results/');
InputDataPath  = fullfile(base_folder,'raw_data/');

cd(fullfile(base_folder))
data_aggregated_already = 1;
if data_aggregated_already
    load(fullfile(InputDataPath,'DATA_2024-05-24.mat'))
else
    DATA = load_raw_in_DATA_per_task(InputDataPath);
end


%%
TGA_results = cell(1,numel(DATA));
%parfor task = 1:numel(TGA_results)
for task = 1:numel(TGA_results)
    data = DATA{task};
    %     data.trial = DATA{task}.trial(1:10);
    %     data.time = DATA{task}.time(1:10);
    %     data.fsample = DATA{task}.fsample;
    %     data.label = DATA{task}.label;

    %% define cfg for TEprepare.m
    cfgTEP = [];

    cfgTEP.toi                 = [min(DATA{task}.time{1,1}),max(DATA{task}.time{1,1})]; % time of interest
    %cfgTEP.sgncmb              = {'A1' 'A2';'A2' 'A1'};  % channels to be analyzed
    cfgTEP.channel             = DATA{task}.label;

    % scanning of interaction delays u
    cfgTEP.predicttimemin_u    = 0;      % minimum u to be scanned
    cfgTEP.predicttimemax_u    = 1000;	  % maximum u to be scanned
    cfgTEP.predicttimestepsize = 200; 	  % time steps between u's to be scanned

    % estimator
    cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

    % ACT estimation and constraints on allowed ACT(autocorelation time)
    cfgTEP.trialselect = 'no';
    cfgTEP.actthrvalue = 10;   % threshold for ACT
    cfgTEP.maxlag      = 10;
    cfgTEP.minnrtrials = 1;   % minimum acceptable number of trials

    % optimizing embedding
    cfgTEP.optimizemethod = 'ragwitz';  % criterion used
    cfgTEP.ragdim         = DATA{task}.dim(1,1):DATA{task}.dim(1,1)*3;       % criterion dimension
    % cfgTEP.ragtaurange    = [0.2 .4];  % range for tau
    % cfgTEP.ragtausteps    = 5;         % steps for ragwitz tau steps
    cfgTEP.ragtaurange    = [.25 .5];
    cfgTEP.ragtausteps    = 5;         % steps for ragwitz tau steps
    cfgTEP.repPred        = 1e3; % size(data.trial{1,1},2)*(1/2);

    % kernel-based TE estimation
    cfgTEP.flagNei = 'Mass';            % neigbour analyse type
    cfgTEP.sizeNei = 4;                 % neigbours to analyse

    % set the level of verbosity of console outputs
    cfgTEP.verbosity = 'info_minor';
    %cfgTEP.verbosity = 'none';

    % define cfg for TE
    cfgTESS = [];

    % use individual dimensions for embedding
    cfgTESS.optdimusage = 'maxdim';
    % cfgTESS.optdimusage = 'indivdim';

    % cfgTESS.dim = DATA{task}.dim(1,1);
    % cfgTESS.tau = DATA{task}.tau(1,1);
    % cfgTESS.dim and cfgTESS.tau

    % statistical and shift testing
    cfgTESS.tail           = 1;
    % cfgTESS.numpermutation = 'findDelay';
    cfgTESS.numpermutation = 1e3;
    % cfgTESS.numpermutation = 1e2;
    cfgTESS.shifttest      = 'no';

    cfgTESS.surrogatetype = 'trialshuffling'; % 'trialreverse','blockreverse2','swapneighbors';
    % cfgTESS.surrogatetype = 'blockreverse2'; % 'blockreverse2','swapneighbors';

    % don't calculate MI additionally to TE
    cfgTESS.MIcalc = 0;

    % results file name
    cfgTESS.fileidout = strcat(OutputDataPath,'out_');

    %% TE compuation, scanning over specified values of u
    TGA_results{task} = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);


    %% Null distribution from surrogate pairs
    % # repetitions of re-analyzing the full data w/ shifted (surrogate) pairs.
    TGA_results{task}.ntrials = size(TGA_results{task}.TEmat,2);
    TGA_results{task}.trial_shuffle_surr_n_per_trial = 3;
    TGA_results{task}.TEmat_sur2 = nan(numel(cfgTEP.channel),0);

    for repeats = 1:TGA_results{task}.trial_shuffle_surr_n_per_trial
        trials2 = randperm(TGA_results{task}.ntrials);
        data2 = data;
        for tr = 1:numel(trials2)
            data2.trial{tr}(2,:) = data.trial{trials2(tr)}(2,:);
        end

        cfgTEP2 = cfgTEP;
        cfgTEP2.predicttimemin_u    = min(TGA_results{task}.cfg.predicttime_u);      % minimum u to be scanned
        cfgTEP2.predicttimemax_u    = max(TGA_results{task}.cfg.predicttime_u);	  % maximum u to be scanned
        if cfgTEP2.predicttimemin_u == cfgTEP2.predicttimemax_u
            cfgTEP2.predicttimemax_u = cfgTEP2.predicttimemin_u + 100;
        end
        cfgTEP2.predicttimestepsize = (cfgTEP2.predicttimemax_u - cfgTEP2.predicttimemin_u); 	  % time steps between u's to be scanned
        cfgTEP2.ragdim = max(TGA_results{task}.cfg.dim);

        TGAtemp = InteractionDelayReconstruction_calculate(cfgTEP2,cfgTESS,data2);
        TGA_results{task}.TGA_surrogate_n{repeats} = TGAtemp;
        TGA_results{task}.TEmat_sur2 = [TGA_results{task}.TEmat_sur2 TGAtemp.TEmat];
        fprintf('%6.2f%% surrogates done.\n',repeats/TGA_results{task}.trial_shuffle_surr_n_per_trial*1e2)
    end
end
save([OutputDataPath 'out_TGA_results_' char(datetime('now','Format','yyyy-MM-dd')) '.mat'],'TGA_results','A','OutputDataPath','InputDataPath');


return

%% Collect in long table for stats.
TEtable = [];
for task=1:numel(TGA_results)
    TEtable = vertcat(TEtable,[DATA{task}.conditions,...
        TGA_results{task}.TEmat' ...
        DATA{task}.dim',...
        DATA{task}.tau'./DATA{task}.fsample,...
        TGA_results{task}.cfg.dim'.*ones(size(TGA_results{task}.TEmat,2),1) ...
        (squeeze(TGA_results{task}.ACT.actvalue(1,:,:))'.*TGA_results{task}.cfg.tau')./DATA{task}.fsample ...
        TGA_results{task}.TEpermvalues(:,1)'.*ones(size(TGA_results{task}.TEmat,2),1) ...
        TGA_results{task}.TEpermvalues(:,6)'.*ones(size(TGA_results{task}.TEmat,2),1) ]);
end
%TGA_results{pp}.cfg.tau.*ones(size(TGA_results{pp}.TEmat,2),1) ...
column_labels = {'Date','Time','Task','Auditory','Visual','epsilon','TEstim','TEpp',...
    'dimstim','dimpp','taustim','taupp','dimhatstim','dimhatpp','tauhatstim','tauhatpp',...
    'permtest12','permtest21','interactdelay12','interactdelay21'};
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
boxplot(TEtable.permtest12,TEtable.Task)
boxplot(TEtable.permtest21,TEtable.Task)
%}


%% Visually check the TE and surrogate distibutions and stats.
for task=1:numel(TGA_results)
    for d=1:size(TGA_results{task}.TEmat_sur,1)
        [c0,n0]=hist(TGA_results{task}.TEmat_sur(d,:));
        [c,n]=hist(TGA_results{task}.TEmat(d,:));
        subplot(2,1,1)
        plot(n0,c0,'--',n,c,'-')
        hold on
    end
    hold off
    legend('stim0->user0','stim->user','user0->stim0','user->stim')
    subplot(2,1,2)
    plot(TGA_results{task}.TEmat')
    fprintf('%8.0f%8.0f%8.0f%8.0f%8.0f%8.3f\n',mean(DATA{task}.conditions))
    fprintf('\n')
    disp(TGA_results{task}.TEpermvalues)
    disp(TGA_results{task}.cfg.dim')
    pause
end

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
    hold off
end
