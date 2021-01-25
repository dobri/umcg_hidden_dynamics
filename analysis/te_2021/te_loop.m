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
    % data
    data = DATA{pp};
    
    %% define cfg for TEprepare.m
    cfgTEP = [];
    
    cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
    %cfgTEP.sgncmb              = {'A1' 'A2';'A2' 'A1'};  % channels to be analyzed
    cfgTEP.channel             = data.label;
    
    % scanning of interaction delays u
    cfgTEP.predicttimemin_u    = 10;      % minimum u to be scanned
    cfgTEP.predicttimemax_u    = 210;	  % maximum u to be scanned
    cfgTEP.predicttimestepsize = 20; 	  % time steps between u's to be scanned
    
    % estimator
    cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)
    
    % ACT estimation and constraints on allowed ACT(autocorelation time)
    cfgTEP.actthrvalue = 100;   % threshold for ACT
    cfgTEP.maxlag      = 100;
    cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials
    
    % optimizing embedding
    cfgTEP.optimizemethod ='ragwitz';  % criterion used
    cfgTEP.ragdim         = 5:9;       % criterion dimension
    %cfgTEP.ragdim         = 3;       % data.dim(1,1) criterion dimension
    cfgTEP.ragtaurange    = [0.75 1.]; % range for tau
    cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
    cfgTEP.repPred        = 100; % round(size(data.trial{1,1},2)*(1/2));      % size(data.trial{1,1},2)*(3/4);
    
    % kernel-based TE estimation
    cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
    cfgTEP.sizeNei = 4;                 % neigbours to analyse
    
    % set the level of verbosity of console outputs
    %cfgTEP.verbosity = 'info_minor';
    cfgTEP.verbosity = 'none';
    
    %% define cfg for TE
    cfgTESS = [];
    
    % use individual dimensions for embedding
    % cfgTESS.optdimusage = 'maxdim';
    cfgTESS.optdimusage = 'indivdim';
    
    %cfgTESS.dim = data.dim(1,1);
    %cfgTESS.tau = data.tau(1,1);
    % cfgTESS.dim and cfgTESS.tau
    
    % statistical and shift testing
    cfgTESS.tail           = 1;
    cfgTESS.numpermutation = 5e3;
    cfgTESS.shifttest      = 'no';
    cfgTESS.alpha          = 1e-3;
    %cfgTESS.shifttesttype  ='TEshift>TE';
    
    cfgTESS.surrogatetype  = 'trialshuffling';
    %cfgTESS.surrogatetype  = 'trialreverse';
    %cfgTESS.surrogatetype  = 'blockreverse2';
    %cfgTESS.surrogatetype  = 'swapneighbors';
    
    % don't calculate MI additionally to TE
    cfgTESS.MIcalc = 0;
    
    % results file name
    cfgTESS.fileidout  = strcat(OutputDataPath,'out_');
    
    %% calculation - scan over specified values for u
    TGA_results{pp} = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);
end
%save([OutputDataPath 'out_TGA_results.mat'],'TGA_results','DATA','OutputDataPath','InputDataPath','TEtable');


return

%% Collect in long table for stats.
TEtable = [];
for pp=1:numel(TGA_results)
    TEtable = vertcat(TEtable,[DATA{pp}.conditions,...
        TGA_results{pp}.TEmat' ...
        DATA{pp}.dim',...
        DATA{pp}.tau'./DATA{pp}.fsample,...
        TGA_results{pp}.cfg.dim'.*ones(size(TGA_results{pp}.TEmat,2),1) ...
        (squeeze(TGA_results{pp}.ACT.actvalue(1,:,:))'.*TGA_results{pp}.cfg.tau')./DATA{pp}.fsample ...
        TGA_results{pp}.TEpermvalues(:,1)'.*ones(size(TGA_results{pp}.TEmat,2),1) ...
        TGA_results{pp}.TEpermvalues(:,6)'.*ones(size(TGA_results{pp}.TEmat,2),1) ]);
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
for pp=1:numel(TGA_results)
    for d=1:size(TGA_results{pp}.TEmat_sur,1)
        [c0,n0]=hist(TGA_results{pp}.TEmat_sur(d,:));
        [c,n]=hist(TGA_results{pp}.TEmat(d,:));
        subplot(2,1,1)
        plot(n0,c0,'--',n,c,'-')
        hold on
    end
    hold off
    legend('stim0->user0','stim->user','user0->stim0','user->stim')
    subplot(2,1,2)
    plot(TGA_results{pp}.TEmat')
    fprintf('%8.0f%8.0f%8.0f%8.0f%8.0f%8.3f\n',mean(DATA{pp}.conditions))
    fprintf('\n')
    disp(TGA_results{pp}.TEpermvalues)
    disp(TGA_results{pp}.cfg.dim')
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
