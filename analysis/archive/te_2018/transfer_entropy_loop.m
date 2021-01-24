% Compute transfer entropy for the single trial and n surrogate trials. Do
% it both ways. The output is the TE values per trial for T->P and P->T
% and the sig of the surrogate test in each case. That's four columns.

% Group comparison? Not necessarily. Just compare the proportion of sig.
% trials. So it's one value for each participant (in training trails). The
% rest of the stats are the same.

%% Setup paths
str=computer;
if strcmp(str,'PCWIN64')
    cd C:\Users\Dobri\Desktop
    addpath('C:\Users\Dobri\Documents\MATLAB\fieldtrip')
    ft_defaults
    addpath('C:\Users\Dobri\Documents\MATLAB\TRENTOOL3\')
    %addpath('C:\Users\Dobri\Desktop\DGD\STAR_data\te')
    %InputDataPath = 'C:\Users\Dobri\Desktop\temp';
    %OutputDataPath = 'groups';
    %raw_data_file='C:\Users\Dobri\Desktop\dgd_larger_data_files\head_bobbing_lv_data_for_spike_sync_surrogates_2017-09-06-16-40.mat'
else
    cd ~/Desktop/c3experiment
    addpath('~/matlab-playwith/fieldtrip/')
    ft_defaults
    addpath('~/matlab-playwith/TRENTOOL3/')
    addpath('~/c3/analysis/te');
    InputDataPath = '~/Desktop/c3experiment/';
    OutputDataPath = fullfile(InputDataPath,'out');
    %raw_data_file='~/Desktop/c3/full_data_set.mat';
end

%%
% Don't search for dimension, hard-code 3, or 2, depending on condition.
% Run it three times (sin pre and post, lorenz pre and post, training).
% So the raw data must be pre-processed into structures per participant and
% stimulus type. Then the output files are saved (just in case) and the
% results re-arranged in a matrix.

fileCell=dir(fullfile(InputDataPath,'Raw*.mat'));
cd(OutputDataPath)
TGA_results_trials=cell(size(fileCell,1),1);

%%
DATA=cell(1,length(fileCell));
for f=1:length(fileCell)
    data=load(fullfile(fileCell(f).folder,fileCell(f).name),'data');
    data=data.data;
    for c=1:size(data.trial,2)
        data.trial{c}(2,:)=smooth(data.trial{c}(2,:),10,'moving');
    end
    DATA{f}=data;
end
clear data

delete(gcp('nocreate'));
poolobj = gcp;
addAttachedFiles(poolobj,{'zscore.m'})
addAttachedFiles(poolobj,{'geornd.m'})
parfor f=1:length(fileCell)
    disp(f)
    % profile on

    % Load
    data=DATA{f};
    data.fsample=mean(data.fsample);
    maxl=Inf;
    for t=1:size(data.time,2)
        maxl=min([maxl size(data.time{t},2)]);
    end
    for t=1:size(data.time,2)
        if maxl==size(data.trial{t},2)
            ind=t;
            break
        end
    end
    for t=1:size(data.time,2)
        data.trial{t}=data.trial{t}(:,1:maxl);
        data.time{t}=data.time{t}(:,1:maxl);
    end
    maxt=max(max(data.time{1}));

    
    %% define cfg for TEprepare.m
    cfgTEP = [];
    cfgTEP.verbosity = 'none';%'info_major';
    
    % data
    cfgTEP.toi                 = [0,maxt]; % time of interest
    cfgTEP.channel             = {'Tutor','PP'};  % channels to be analyzed
    
    % scanning of interaction delays u
    cfgTEP.predicttimemin_u    = 20;      % minimum u to be scanned
    cfgTEP.predicttimemax_u    = 520;	  % maximum u to be scanned
    cfgTEP.predicttimestepsize = 50; 	  % time steps between u's to be scanned
    
    % estimator
    cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)
    
    % ACT estimation and constraints on allowed ACT(autocorelation time)     <-
    % The faughtdoodles keep switching unit convention between ms and sample.
    cfgTEP.actthrvalue = 100;   % threshold for ACT
    cfgTEP.maxlag      = 100;
    cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials
    cfgTEP.trialselect = 'no';
    
    % optimizing embedding                                                   <-
    cfgTEP.optimizemethod ='ragwitz';  % criterion used
    cfgTEP.ragdim         = data.dim{1}(1);       % criterion dimension
    cfgTEP.ragtaurange    = [0.2 0.5]; % range for tau
    cfgTEP.ragtausteps    = 5;         % steps for ragwitz tau steps
    cfgTEP.repPred        = round(size(data.trial{1},2)*(1/2));    % size(data.trial{1,1},2)*(3/4);   %?
    
    % kernel-based TE estimation
    cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
    cfgTEP.sizeNei = 4;                 % neigbours to analyse
    
    cfgTEP.ensemblemethod = 'no';
    cfgTEP.outputpath = OutputDataPath;
    
    cfgTEP.TEparallel.parON   = 'no'; %'yes'; % for paralell computing tool.
    cfgTEP.TEparallel.workers =   16; % number of workers for parallel computing.
    
    
    
    %% define cfg for TEsurrogatestats_ensemble.m
    
    cfgTESS = [];
    
    % use individual dimensions for embedding
    cfgTESS.optdimusage = 'indivdim';
    
    % statistical and shift testing
    cfgTESS.tail           = 1;
    cfgTESS.numpermutation = 5e2;
    cfgTESS.tail           = 1;
    cfgTESS.surrogatetype  = 'blockreverse4';
    % cfgTESS.surrogatetype  = 'trialshuffling';
    % cfgTESS.surrogatetype  = 'blockresampling'; %X % trialreverse
    % cfgTESS.surrogatetype  = 'blockreverse1'; % V
    % cfgTESS.surrogatetype  = 'blockreverse2';
    
    shuffle_n = 1e3;
    
    cfgTESS.shifttest      = 'no';
    cfgTESS.shifttesttype  = 'TEshift>TE';
    
    
    % results file name
    cfgTESS.fileidout  = fullfile(OutputDataPath,fileCell(f).name);
    
    %% Shoot
    TGA_results = basic_TEpermtest(cfgTEP,cfgTESS,data,strcat(cfgTEP.outputpath,fileCell(f).name(1:end-4)),shuffle_n);
    TGA_results_trials{f} = TGA_results;
    % for tt=1:8;surtest=sum(TGA_results.TEmat_shuffle(:,tt)<TGA_results.TEmat(tt))/size(TGA_results.TEmat_shuffle,1)*100;disp(surtest);end
end
