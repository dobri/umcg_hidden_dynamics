%% set paths

base_folder = '/home/lilabuntu/Desktop/c3/umcg/';
addpath('~/Desktop/TRENTOOL3')
addpath('~/Desktop/TRENTOOL3/fieldtrip');
ft_defaults;

%% define data paths

OutputDataPath = fullfile(base_folder,'te_2021/results/');
InputDataPath  = fullfile(base_folder,'data/training_processed/training');

%load(fullfile(InputDataPath,'lorenz_1-2_45ms.mat'));
%DATA = convert_csv_to_ft_and_te(InputDataPath);
cd(fullfile(base_folder,'te_2021'))

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
    cfgTEP.predicttimemin_u    = 40;      % minimum u to be scanned
    cfgTEP.predicttimemax_u    = 160;	  % maximum u to be scanned
    cfgTEP.predicttimestepsize = 20; 	  % time steps between u's to be scanned
    
    % estimator
    cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)
    
    % ACT estimation and constraints on allowed ACT(autocorelation time)
    cfgTEP.actthrvalue = 100;   % threshold for ACT
    cfgTEP.maxlag      = 100;
    cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials
    
    % optimizing embedding
    cfgTEP.optimizemethod ='ragwitz';  % criterion used
    cfgTEP.ragdim         = 2:9;       % criterion dimension
    %cfgTEP.ragdim         = 3;       % data.dim(1,1) criterion dimension
    cfgTEP.ragtaurange    = [0.1 0.5]; % range for tau
    cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
    cfgTEP.repPred        = round(size(data.trial{1,1},2)*(1/2));      % size(data.trial{1,1},2)*(3/4);
    
    % kernel-based TE estimation
    cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
    cfgTEP.sizeNei = 4;                 % neigbours to analyse
    
    % set the level of verbosity of console outputs
    %cfgTEP.verbosity = 'info_minor';
    cfgTEP.verbosity = 'none';
    
    %% define cfg    
    cfgTESS = [];
    
    % use individual dimensions for embedding
    cfgTESS.optdimusage = 'maxdim';%'indivdim';
    
    %cfgTESS.dim = data.dim(1,1);
    %cfgTESS.tau = data.tau(1,1);
    
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
%save([OutputDataPath 'out_TGA_results.mat'],'TGA_results','DATA','OutputDataPath','InputDataPath','Surr');


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


Surr = zeros(40,2,3);
Surrcount= zeros(40,1,3);
for f=1:numel(DATA)
    switch DATA{f}.conditions(1,3)
        case 10
            task = 1;
        case 25
            task = 2;
        case 30
            task = 3;
    end
    Surrcount(1:numel(TGA_results{f}.TEmat(1,:)'),1,task) = Surrcount(1:numel(TGA_results{f}.TEmat(1,:)'),1,task)+1;
    Surr(1:numel(TGA_results{f}.TEmat(1,:)'),1,task) = Surr(1:numel(TGA_results{f}.TEmat(1,:)'),1,task)+TGA_results{f}.TEmat(1,:)';
    Surr(1:numel(TGA_results{f}.TEmat(1,:)'),2,task) = Surr(1:numel(TGA_results{f}.TEmat(1,:)'),2,task)+TGA_results{f}.TEmat(2,:)';
end
for task=1:3
    Surr(:,1,task) = Surr(:,1,task)./Surrcount(:,1,task);
    Surr(:,2,task) = Surr(:,2,task)./Surrcount(:,1,task);
    Surr(:,3,task) = task;
end
Surr = reshape(permute(Surr,[1 3 2]),[],3);
for task=1:3
    subplot(1,3,task)
    plot(Surr(Surr(:,3)==task,1),'-ob')
    hold on
    plot(Surr(Surr(:,3)==task,2),'-or')
    hold off
end


return

%% optional: perform a post hoc correction for cascade effects and simple common drive effects

cfgGA = [];

cfgGA.threshold = 3;
cfgGA.cmc       = 1;

TGA_results_GA = TEgraphanalysis(cfgGA,TGA_results);

%save([OutputDataPath 'Lorenz_1->2_TGA_results_analyzed_GA.mat'],'TGA_results_GA');


%% plotting

load(fullfile(InputDataPath,'lorenz_layout.mat'))


cfgPLOT = [];

cfgPLOT.layout        = lay_Lorenz; 		% see fieldtrip's ft_prepare_layout.m
cfgPLOT.electrodes    = 'highlights';
cfgPLOT.statstype     = 1;   		% 1: corrected; 2:uncorrected; 3: 1-pval; 4:rawdistance
cfgPLOT.alpha         = 0.05;
cfgPLOT.arrowpos      = 1;
cfgPLOT.showlabels    = 'yes';
cfgPLOT.electrodes    = 'on';
cfgPLOT.hlmarker      = 'o';
cfgPLOT.hlcolor       = [0 0 0];
cfgPLOT.hlmarkersize  = 4;
cfgPLOT.arrowcolorpos = [1 0 0];

figure;
TEplot2D(cfgPLOT,TGA_results_GA)
