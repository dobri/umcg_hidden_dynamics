%% set paths

base_folder = '~/Desktop/TRENTOOL3/working_dir';
cd(base_folder)
addpath('~/Desktop/TRENTOOL3')
addpath('~/Desktop/TRENTOOL3/fieldtrip');
ft_defaults;

%% define data paths

OutputDataPath = fullfile(base_folder,'results/');
InputDataPath  = fullfile(base_folder,'../../TRENTOOL3_exampledata/Lorenz_2_systems/');

%load(fullfile(InputDataPath,'lorenz_1-2_45ms.mat'));
DATA = convert_csv_to_ft_and_te('/home/lilabuntu/Desktop/c3/umcg/data/training_processed/training');
data = DATA{1};

%% define cfg for TEprepare.m

cfgTEP = [];


% data
cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
%cfgTEP.sgncmb              = {'A1' 'A2';'A2' 'A1'};  % channels to be analyzed
cfgTEP.channel             = {'A1','A2'}; 

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 44;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 46;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 1; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = 1000;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:9;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100;      % size(data.trial{1,1},2)*(3/4);

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse

% set the level of verbosity of console outputs
%cfgTEP.verbosity = 'info_minor';
cfgTEP.verbosity = 'none';

%% define cfg for TEsurrogatestats_ensemble.m

cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 5e2;
cfgTESS.shifttest      = 'no';
cfgTESS.shifttesttype  ='TEshift>TE';
cfgTESS.surrogatetype  = 'trialshuffling';
%cfgTESS.surrogatetype  = 'trialreverse';
%cfgTESS.surrogatetype  = 'blockreverse2';
%cfgTESS.surrogatetype  = 'swapneighbors';

% don't calculate MI additionally to TE
cfgTESS.MIcalc = 0;

% results file name
cfgTESS.fileidout  = strcat(OutputDataPath,'Lorenzdata_1->2_');

%% calculation - scan over specified values for u

TGA_results = cell(4,1);
parfor r = 1:numel(TGA_results)
    TGA_results{r} = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);
end
%save([OutputDataPath 'Lorenz_1->2_TGA_results.mat'],'TGA_results');

for p=1:numel(TGA_results)
    for r=1:2
        [c0,n0]=hist(TGA_results{p}.TEmat_sur(r,:));
        [c,n]=hist(TGA_results{p}.TEmat(r,:));
        subplot(2,1,r);plot(n0,c0,'-b',n,c,'--r')
    end
    pause
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
